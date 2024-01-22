#pragma once

#include <chrono>
#include <vector>
#include <iostream>
#include <thread>
#include <fstream>


//#define PROFILING 1

#if PROFILING

#define PROFILE_SCOPE(name) ::benchmark::Timer timer##__LINE__(name);
#define PROFILE_FUNC() PROFILE_SCOPE(BOOST_CURRENT_FUNCTION)

#else

#define PROFILE_SCOPE(name)
#define PROFILE_FUNC()

#endif



namespace benchmark {

	struct BenchmarkData {
		long long start_time{ 0 };
		long long exec_time{ 0 };
		const char *func_name{ "" };
		uint32_t thread_id{ 0 };
	};


	struct global_timer {


		static global_timer &get_inst() {
			static global_timer instance{};
			return instance;
		}

		static void push(BenchmarkData data) {
			get_inst().benchmarks.push_back(data);
		}

		//* viewable in about:tracing *//

		static void write_to_file(const std::string &name) {
			auto &inst = get_inst();
			std::ofstream ostr(name);

			ostr << "{\"otherData\": {},\"traceEvents\":[";

			bool first = true;

			for (auto &b : inst.benchmarks) {

				if (first) {
					first = false;
				}
				else {
					ostr << ",";
				}

				std::string name = b.func_name;
				std::replace(name.begin(), name.end(), '"', '\'');

				ostr << "{";
				ostr << "\"cat\":\"function\",";
				ostr << "\"dur\":" << b.exec_time << ',';
				ostr << "\"name\":\"" << name << "\",";
				ostr << "\"ph\":\"X\",";
				ostr << "\"pid\":0,";
				ostr << "\"tid\":" << b.thread_id << ",";
				ostr << "\"ts\":" << b.start_time;

				ostr << "}";
			}

			ostr << "]}";
		}

		std::vector<BenchmarkData> benchmarks{};

		global_timer(global_timer const &) = delete;
		void operator=(global_timer const &) = delete;


	private:

		global_timer() = default;

	};


	struct Timer {

		Timer(const char *name) {
			start(name);
		}

		~Timer() {
			stop();
		}

		inline void start(const char *func_name = nullptr) {
			start_time = std::chrono::high_resolution_clock::now();
			name = func_name;
		}

		inline void stop() {
			auto end_time = std::chrono::high_resolution_clock::now();

			auto start = std::chrono::time_point_cast<std::chrono::microseconds>(start_time).time_since_epoch().count();
			auto end = std::chrono::time_point_cast<std::chrono::microseconds>(end_time).time_since_epoch().count();

			uint32_t thread_id = (uint32_t)std::hash<std::thread::id>{}(std::this_thread::get_id());

			auto data = BenchmarkData{
				.start_time = start,
					.exec_time = end - start,
					.func_name = name,
					.thread_id = thread_id,
			};

			global_timer::push(std::move(data));
		}

		const char *name;
		std::chrono::time_point<std::chrono::high_resolution_clock> start_time{};
	};


}
