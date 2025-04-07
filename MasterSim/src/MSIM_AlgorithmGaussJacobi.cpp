#include "MSIM_AlgorithmGaussJacobi.h"

#include <IBK_assert.h>

#include "MSIM_MasterSim.h"

// hrkim 
// #include <omp.h>
#include <thread>
// #include <vector>

#include <IBK_messages.h>


namespace MASTER_SIM {

void AlgorithmGaussJacobi::init() {
	const char * const FUNC_ID = "[AlgorithmGaussJacobi::init]";
	// check for valid parameters
	if (m_master->m_project.m_maxIterations > 1)
		throw IBK::Exception("Gauss-Jacobi-Algorithm is always used without iteration. The maxIterations parameter must be set to 1.", FUNC_ID);
}


AbstractAlgorithm::Result AlgorithmGaussJacobi::doStep() {
	const char * const FUNC_ID = "[AlgorithmGaussJacobi::doStep]";

	// to make things simpler, let's just use fmi2status variables
	IBK_STATIC_ASSERT((int)fmiOK == (int)fmi2OK);

	// master and FMUs are expected to be at current time point t

	// all slave output variables are expected to be in sync with internal states of slaves
	// i.e. cacheOutputs() has been called successfully on all slaves

	// global variable array is expected to be in sync with all slaves

	// ** algorithm start **

	// loop over all cycles
		// #pragma omp parallel for //reduction(max:maxSlaveTime)
		for (unsigned int c=0; c<m_master->m_cycles.size(); ++c) {
			const MasterSim::Cycle & cycle = m_master->m_cycles[c];
			
			++m_nIterations; // add one iteration
			
			
			// hrkim
			// double maxSlaveTime = 0.0;
			// double parallel_start = omp_get_wtime();
			
			// hrkim end
			// loop over all slaves
			// #pragma omp parallel for //schedule(dynamic,1)
			
			// hrkim : std::vector to hold threads for this cycle
			// std::vector<std::thread> threads;

			for (unsigned int s=0; s<cycle.m_slaves.size(); ++s) {
				// hrkim
				// double slave_start = omp_get_wtime();
				// threads.push_back(std::thread([&, s]() {


				AbstractSlave * slave = cycle.m_slaves[s];
				
				
				// update input variables in all slaves, using variables from time t
				m_master->updateSlaveInputs(slave, m_master->m_realyt, m_master->m_intyt, m_master->m_boolyt, m_master->m_stringyt, false);

				
				// hrkim
				double elapsed_time;
				// advance slave
				m_timer.start();
				// double slave_start = omp_get_wtime();
				int res = slave->doStep(m_master->m_h, true);

				// double slave_time = omp_get_wtime() - slave_start;
				// hrkim
				elapsed_time = 1e-3*m_timer.stop();
				
				
				m_master->m_statSlaveEvalTimes[slave->m_slaveIndex] += elapsed_time; // add elapsed time in seconds
				++m_master->m_statSlaveEvalCounters[slave->m_slaveIndex];
				
				if (res != fmi2OK)
				throw IBK::Exception(IBK::FormatString("Error in doStep() call of FMU slave '%1'").arg(slave->m_name), FUNC_ID);
				
				// slave is now at time level t + h and its outputs are updated accordingly
				// sync results into vector with newly computed quantities
				m_master->syncSlaveOutputs(slave, m_master->m_realytNext, m_master->m_intytNext, m_master->m_boolytNext, m_master->m_stringytNext, false);
				
				// hrkim
				// double slave_time = omp_get_wtime() - slave_start;
				// maxSlaveTime = (slave_time > maxSlaveTime) ? slave_time : maxSlaveTime;
				// })); // hrkim : end of thread
			}
			// hrkim : thread join : wait until the end of all slave's doStep()
			// for (auto & t : threads) {
			// 	if (t.joinable())
			// 		t.join();
			// }
			// hrkim 
			// double parallel_end = omp_get_wtime();
			// double overallParallelTime = parallel_end - parallel_start;
			
			// double overheadTime = overallParallelTime - maxSlaveTime;
			
			// accumulatedOverhead += overheadTime;
			// accumulatedParallelTime += overallParallelTime;
			
			// if (cnt == 9997)
			// {ECU/script/EXP2$ ./exp2_ECU3.sh

				// std::cout << cnt << " : 전체 accumulated 병렬 실행 시간 = " << accumulatedParallelTime 
				// << " sec, 최대 doStep() 실행 시간 = " << maxSlaveTime 
				// << " sec, 추정 오버헤드 = " << overheadTime << " sec\n"
				// << "누적 오버헤드 = " << accumulatedOverhead << " sec\n";
				// }
				// cnt++;
				// hrkim end
			// }
			
			// ** algorithm end **
			
			// m_XXXyt     -> still values at time point t
			// m_XXXytNext -> values at time point t + h
		}
		// hrkim : pragma end brace
		// }
			return R_CONVERGED; // no other option since we don't iterate
		}
		

} // namespace MASTER_SIM
