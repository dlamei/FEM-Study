#include "FEM_ocl.h"

#include "utils.h"

#include <CL/cl.hpp>
#include <vector>


#pragma comment(lib, "OpenCL.lib")



namespace cl_utils {
	struct Context {
		cl::Platform platform{};
		cl::Device device{};
		cl::Context context{};

		cl::CommandQueue main_queue{};

		static Context init() {
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

			cl::Context context({ device });
			cl::CommandQueue main_queue(context, device);

			return Context{
				.platform = platform,
				.device = device,
				.context = context,
				.main_queue = main_queue,
			};
		}
	};

	enum class MemUsage {
		READ_WRITE = CL_MEM_READ_WRITE,
		READ = CL_MEM_READ_ONLY,
		WRITE = CL_MEM_WRITE_ONLY,
	};

	struct Buffer {
		usize size;
		cl::Buffer data;

		static Buffer empty(const Context &c, usize size, MemUsage usage = MemUsage::READ_WRITE) {
			cl::Buffer data(c.context, (u64)MemUsage::READ_WRITE, size);

			return Buffer{
				.size = size,
				.data = data,
			};
		}

		static Buffer init(const Context &c, void *data_ptr, usize size, MemUsage usage = MemUsage::READ_WRITE) {
			auto b = Buffer::empty(c, size, usage);
			c.main_queue.enqueueWriteBuffer(b.data, CL_TRUE, 0, size, data_ptr);
			return b;
		}

		void to_cpu(const Context &c, void *storage_ptr, usize size = -1) {
			if (size == -1) size = this->size;
			assert(this->size >= size);
			c.main_queue.enqueueReadBuffer(data, CL_TRUE, 0, size, storage_ptr);
		}
	};

}

namespace fem_ocl {


	void solve()
	{
		auto c = cl_utils::Context::init();

		cl::Program::Sources sources;

		std::string kernel_code = R"(
		void kernel simple_add(global const int *A, global const int *B, global int *C) {
			int id = get_global_id(0);
			C[id] = A[id] + B[id];
		}
		)";


		sources.push_back({ kernel_code.c_str(), kernel_code.length() });


		cl::Program program(c.context, sources);
		if (program.build({ c.device }) != CL_SUCCESS) {
			auto msg = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(c.device);
			assert(false, msg.c_str());
		}


		i32 a[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
		i32 b[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

		auto b1 = cl_utils::Buffer::init(c, a, 10 * sizeof(i32));
		auto b2 = cl_utils::Buffer::init(c, b, 10 * sizeof(i32));
		auto b3 = cl_utils::Buffer::empty(c, 10 * sizeof(i32));


		cl::make_kernel<cl::Buffer, cl::Buffer, cl::Buffer> simple_add(cl::Kernel(program, "simple_add"));
		cl::EnqueueArgs eargs(c.main_queue, cl::NullRange, cl::NDRange(10), cl::NullRange);
		simple_add(eargs, b1.data, b2.data, b3.data).wait();

		i32 res[10] = {};
		b3.to_cpu(c, res);

		printf("result: \n");

		for (i32 i = 0; i < 10; i++) {
			printf("%i ", res[i]);
		}
		printf("\n");

	}


}
