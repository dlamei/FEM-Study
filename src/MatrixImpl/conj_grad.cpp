#include "conj_grad.h"

#include "sparse_matrix.h"
#include "../utils.h"

#include <CL/cl.hpp>
#include <vector>


#pragma comment(lib, "OpenCL.lib")

namespace cl_utils {

	struct State {
		cl::Platform platform{};
		cl::Device device{};
		cl::Context context{};

		cl::CommandQueue queue{};

		static State init() {
			std::vector<cl::Platform> platforms;
			cl::Platform::get(&platforms);
			assert(platforms.size(), "no platform found");
			cl::Platform platform = platforms[0];

			std::vector<cl::Device> devices;
			platform.getDevices(CL_DEVICE_TYPE_ALL, &devices);
			assert(devices.size(), "no device found!");
			cl::Device device = devices[0];

			std::cout << "using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << "\n";
			std::cout << "using device: " << device.getInfo<CL_DEVICE_NAME>() << "\n";
			std::cout << "max workgroup size: " << CL_DEVICE_MAX_WORK_GROUP_SIZE << "\n";

			cl::Context context({ device });
			cl::CommandQueue queue(context, device);

			return State{
				.platform = platform,
				.device = device,
				.context = context,
				.queue = queue,
			};
		}
	};

	struct SparseMatrix {
		index_t n_rows{ 0 }, n_cols{ 0 };
		cl::Buffer vals{};
		cl::Buffer col_indx{};
		cl::Buffer row_ptr{};

		static SparseMatrix read_only_from_cpu(const State &s, sparse_matrix::Matrix &m) {
			index_t n_rows = m.rows;
			index_t n_cols = m.cols;

			auto mem_flag = CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR;
			auto vals = cl::Buffer(s.context, mem_flag, m.vals.size() * sizeof(scalar), m.vals.data());
			auto col_indx = cl::Buffer(s.context, mem_flag, m.col_indx.size() * sizeof(index_t), m.col_indx.data());
			auto row_ptr = cl::Buffer(s.context, mem_flag, m.row_ptr.size() * sizeof(index_t), m.row_ptr.data());

			return SparseMatrix{
				.n_rows = n_rows,
				.n_cols = n_cols,
				.vals = vals,
				.col_indx = col_indx,
				.row_ptr = row_ptr,
			};
		}
	};

	struct Vector {
		usize size{};
		cl::Buffer buf{};

		static Vector empty(const State &s, usize size) {
			auto buf = cl::Buffer(s.context, CL_MEM_READ_WRITE, size);

			return Vector{
				.size = size,
				.buf = buf,
			};
		}

		static Vector from_ptr(const State &s, void *data, usize size) {
			auto buf = cl::Buffer(s.context, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, size, data);

			return Vector{
				.size = size,
				.buf = buf,
			};
		}

		template<typename T>
		std::vector<T> to_cpu(const State &s) {
			usize count = size / sizeof(T);
			assert(count);

			std::vector<T> data(count, 0);
			s.queue.enqueueReadBuffer(buf, CL_TRUE, 0, size, data.data());
			return data;
		}

	};

	inline cl::EnqueueArgs enqueue(cl::CommandQueue &queue, index_t glob_size) {
		return cl::EnqueueArgs(queue, cl::NullRange, cl::NDRange(glob_size), cl::NullRange);
	}
}


namespace fem_ocl {

	const std::string KERNEL_SOURCE = {
		#include "conj_grad.ocl" 
	};

	typedef cl::make_kernel<index_t, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer> sqr_mat_vec_mul_kernel;
	typedef cl::make_kernel<scalar, cl::Buffer, cl::Buffer, cl::Buffer> vec_add_mul_kernel;
	typedef cl::make_kernel<scalar, cl::Buffer, scalar, cl::Buffer> vec_add_mul_assign_kernel;
	typedef cl::make_kernel<index_t, cl::Buffer, cl::Buffer> vec_copy_kernel;
	typedef cl::make_kernel<index_t, scalar, cl::Buffer> vec_fill_kernel;
	typedef cl::make_kernel<index_t, scalar, scalar, cl::Buffer, cl::Buffer, cl::Buffer, cl::Buffer> conj_grad_step_kernel;

	struct ConjGradKernels {
		cl::Program program;

		sqr_mat_vec_mul_kernel sqr_mat_vec_mul;
		vec_add_mul_kernel vec_add_mul;
		vec_add_mul_assign_kernel vec_add_mul_assign;
		vec_copy_kernel vec_copy;
		vec_fill_kernel vec_fill;
		conj_grad_step_kernel conj_grad_step;

		static ConjGradKernels load_kernels(const cl_utils::State &state) {

			cl::Program::Sources sources;
			sources.push_back({ KERNEL_SOURCE.c_str(), KERNEL_SOURCE.length() });
			cl::Program program(state.context, sources);

			if (program.build({ state.device }) != CL_SUCCESS) {
				auto msg = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(state.device);
				assert(false, msg.c_str());
			}


			auto sqr_mat_vec_mul = sqr_mat_vec_mul_kernel(cl::Kernel(program, "sqr_mat_vec_mul"));
			auto vec_add_mul = vec_add_mul_kernel(cl::Kernel(program, "vec_add_mul"));
			auto vec_add_mul_assign = vec_add_mul_assign_kernel(cl::Kernel(program, "vec_add_mul_assign"));
			auto vec_copy = vec_copy_kernel(cl::Kernel(program, "vec_copy"));
			auto vec_fill = vec_fill_kernel(cl::Kernel(program, "vec_fill"));
			auto conj_grad_step = conj_grad_step_kernel(cl::Kernel(program, "conj_grad_step"));


			return ConjGradKernels{
				.program = program,
				.sqr_mat_vec_mul = sqr_mat_vec_mul,
				.vec_add_mul = vec_add_mul,
				.vec_add_mul_assign = vec_add_mul_assign,
				.vec_copy = vec_copy,
				.vec_fill = vec_fill,
				.conj_grad_step = conj_grad_step,
			};
		}
	};

	void print(const std::vector<scalar> &vec) {
		for (const auto s : vec) {
			printf("%f, ", s);
		}
		printf("\n");
	}

	scalar dot(const std::vector<scalar> &v1, const std::vector<scalar> &v2) {
		db_assert(v1.size() == v2.size());

		scalar sum = 0;
		for (isize i = 0; i < v1.size(); i++) {
			sum += v1[i] * v2[i];
		}

		return sum;
	}

	// iter solve Ax = b
	std::vector<scalar> run_conj_grad(sparse_matrix::Matrix &_A, std::vector<scalar> &_b) {
		auto s = cl_utils::State::init();
		auto cg = ConjGradKernels::load_kernels(s);

		const index_t n_rows = _A.rows;
		const auto eargs = cl_utils::enqueue(s.queue, n_rows);

		auto tmp = cl_utils::Vector::empty(s, n_rows * sizeof(scalar));

		auto A = cl_utils::SparseMatrix::read_only_from_cpu(s, _A);
		auto b = cl_utils::Vector::from_ptr(s, _b.data(), n_rows * sizeof(scalar));

		//init x = b
		auto x = cl_utils::Vector::from_ptr(s, _b.data(), n_rows * sizeof(scalar));

		// r_0 = b - Ax
		auto r = cl_utils::Vector::empty(s, n_rows * sizeof(scalar));
		cg.sqr_mat_vec_mul(eargs, n_rows, A.vals, A.col_indx, A.row_ptr, x.buf, tmp.buf).wait(); // tmp = Ax
		cg.vec_add_mul(eargs, -1, b.buf, tmp.buf, r.buf).wait(); // r = b - Ax

		// p_0 = r_0
		auto p = cl_utils::Vector::empty(s, n_rows * sizeof(scalar));
		cg.vec_copy(eargs, n_rows, r.buf, p.buf).wait();

		// pre alloc
		auto Ap = cl_utils::Vector::empty(s, n_rows * sizeof(scalar));

		std::vector<scalar> cpu_r(n_rows, 0);
		std::vector<scalar> cpu_p(n_rows, 0);
		std::vector<scalar> cpu_Ap(n_rows, 0);

		s.queue.enqueueReadBuffer(r.buf, CL_TRUE, 0, n_rows * sizeof(scalar), cpu_r.data());
		scalar old_rTr = dot(cpu_r, cpu_r);

		for (u32 i = 0; i < 10000; i++) {

			cl_event mul_ev = cg.sqr_mat_vec_mul(eargs, n_rows, A.vals, A.col_indx, A.row_ptr, p.buf, Ap.buf)();

			s.queue.enqueueReadBuffer(Ap.buf, CL_TRUE, 0, n_rows * sizeof(scalar), cpu_Ap.data());
			s.queue.enqueueReadBuffer(p.buf, CL_TRUE, 0, n_rows * sizeof(scalar), cpu_p.data());

			scalar pTAp = dot(cpu_p, cpu_Ap);
			scalar alpha = old_rTr / pTAp;

			// x = x + alpha * p
			cg.vec_add_mul_assign(eargs, 1, x.buf, alpha, p.buf);
			cg.vec_add_mul_assign(eargs, 1, r.buf, -alpha, Ap.buf);

			s.queue.enqueueReadBuffer(r.buf, CL_TRUE, 0, n_rows * sizeof(scalar), cpu_r.data());

			scalar new_rTr = dot(cpu_r, cpu_r);
			scalar beta = new_rTr / old_rTr;

			// p = r + beta * p
			cg.vec_add_mul_assign(eargs, beta, p.buf, 1, r.buf);

			s.queue.finish();

			old_rTr = new_rTr;

			if (old_rTr < 0.1) {
				printf("equilibrium reached!");
				break;
			}

			if (i % 100 == 0) {
				printf("residium norm: %f\n", std::sqrt(old_rTr));
			}
		}


		return x.to_cpu<scalar>(s);
	}


	//void conj_grad() {

	//	std::vector<Triplet> triplets{};
	//	triplets.push_back({ 0, 0, 4 });
	//	triplets.push_back({ 0, 1, 12 });
	//	triplets.push_back({ 0, 2, -16 });
	//	triplets.push_back({ 1, 0, 12 });
	//	triplets.push_back({ 1, 1, 37 });
	//	triplets.push_back({ 1, 2, -43 });
	//	triplets.push_back({ 2, 0, -16 });
	//	triplets.push_back({ 2, 1, -43 });
	//	triplets.push_back({ 2, 2, 98 });

	//	auto m = sparse_matrix::Matrix::from_triplets(3, 3, &triplets);

	//	std::vector<scalar> b = { 1, 1, 1 };

	//	run_conj_grad(m, b);
	//}

}



TEST(ocl_sqr_mat_vec_mul, {

		std::vector<Triplet> triplets{};
		triplets.push_back(Triplet{.row = 0, .col = 0, .val = 1 });
		triplets.push_back(Triplet{ 0, 1, 2 });
		triplets.push_back(Triplet{ 0, 2, 3 });
		triplets.push_back(Triplet{ 1, 0, 4 });
		triplets.push_back(Triplet{ 1, 1, 5 });
		triplets.push_back(Triplet{ 1, 2, 6 });
		triplets.push_back(Triplet{ 2, 0, 7 });
		triplets.push_back(Triplet{ 2, 1, 8 });
		triplets.push_back(Triplet{ 2, 2, 9 });

		std::vector<scalar> v1 = { 10, 11, 12 };
		std::vector<scalar> v2 = { 1, 2, 3 };
		auto sparse = sparse_matrix::Matrix::from_triplets(3, 3, &triplets);

		index_t n_rows = sparse.rows;


		auto s = cl_utils::State::init();
		auto cg = fem_ocl::ConjGradKernels::load_kernels(s);


		auto cl_mat = cl_utils::SparseMatrix::read_only_from_cpu(s, sparse);
		auto cl_v1 = cl_utils::Vector::from_ptr(s, v1.data(), v1.size() * sizeof(scalar));
		auto cl_v2 = cl_utils::Vector::from_ptr(s, v2.data(), v2.size() * sizeof(scalar));

		auto out1 = cl_utils::Vector::empty(s, v1.size() * sizeof(scalar));
		auto out2 = cl_utils::Vector::empty(s, v1.size() * sizeof(scalar));


		auto qargs = cl_utils::enqueue(s.queue, n_rows);

		// out1 = A * v1
		cg.sqr_mat_vec_mul(qargs, n_rows, cl_mat.vals, cl_mat.col_indx, cl_mat.row_ptr, cl_v1.buf, out1.buf);
		// out2 = v1 - v2
		cg.vec_sub(qargs, n_rows, cl_v1.buf, cl_v2.buf, out2.buf);

		s.queue.finish();

		auto res1 = out1.to_cpu<scalar>(s);
		auto res2 = out2.to_cpu<scalar>(s);

		test_assert(cmp_scalar(68, res1[0]));
		test_assert(cmp_scalar(167, res1[1]));
		test_assert(cmp_scalar(266, res1[2]));

		test_assert(cmp_scalar(9, res2[0]));
		test_assert(cmp_scalar(9, res2[1]));
		test_assert(cmp_scalar(9, res2[2]));

		return {};
	});
