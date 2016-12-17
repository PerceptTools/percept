// Copyright 2014 Sandia Corporation. Under the terms of
// Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef percept_RunEnvironment_hpp
#define percept_RunEnvironment_hpp

/// copied and edited from stk_util/use_cases/UseCaseEnvironment


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/BroadcastArg.hpp>
#include <stk_util/util/Writer_fwd.hpp>
#include <stk_util/diag/Timer.hpp>
#include <stk_util/environment/OutputLog.hpp>
#include <stk_util/environment/ReportHandler.hpp>
#include <stk_util/environment/RuntimeWarning.hpp>
#include <stk_util/environment/RuntimeDoomed.hpp>

#include <iosfwd>

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

  namespace percept {


    enum LogMask {
      LOG_ALWAYS          = stk::LOG_ALWAYS,
      LOG_TRACE           = stk::LOG_TRACE,
      LOG_TRACE_STATS     = stk::LOG_TRACE_STATS,
      LOG_TRACE_SUB_CALLS = stk::LOG_TRACE_SUB_CALLS,
      LOG_MEMBERS         = stk::LOG_MEMBERS,

      LOG_APPLICATION     = 0x0000010  // use this as the base for additional masks
      //LOG_SEARCH          = 0x0000010,
      //LOG_TRANSFER        = 0x0000020,
      //LOG_TIMER           = 0x0000040

    };

    /**
     * @brief Class <b>message_type</b> ...
     *
     */
    enum message_type {
      MSG_WARNING = stk::MSG_WARNING,
      MSG_FATAL   = stk::MSG_DOOMED,
      MSG_INFORMATION,
      MSG_EXCEPTION,
      MSG_PARALLEL_EXCEPTION
    };


    /**
     * @brief Class <b>type</b> ...
     *
     */
    enum message_throttle_type {
      MSG_APPLICATION = stk::MSG_APPLICATION,
      MSG_TIME_STEP
    };

    enum TimerSetMask {
      TIMER_MESH     = 0x00000001,		///< Enable mesh timers
      //      TIMER_MESH_IO  = 0x00000002,		///< Enable mesh I/O timers
      //      TIMER_SEARCH   = 0x00000004,		///< Enable search timers
      //    TIMER_TRANSFER = 0x00000008,		///< Enable transfer timers
      TIMER_ALL      = 0xFFFFFFFF,		///< Force timer to be active

      TIMER_FORCE    = 0x00000000		///< Force timer to be active
    };

    std::ostream &out();                ///< Normal output stream
    std::ostream &dout();               ///< Diagnostic output stream
    std::ostream &pout();               ///< Per-processor output stream (See RuntimeDeferredx)
    std::ostream &tout();               ///< Regression test textual output stream

    std::ostream &dwout();              ///< Diagnostic writer stream

    // dw() definition
    stk::diag::Writer &dw();
#define DWENDL stk::diag::dendl

    stk::diag::TimerSet &timerSet();

    stk::diag::Timer &timer();

    void my_report_handler(const char *message, int type);

    // this little class is simply here to force an ordering in the ~RunEnvironment() dtor, so that
    //   parallel_machine_finalize gets called after m_comm is destroyed; but, it still only works
    //   if this is invoked after returning from an enclosing block of RunEnvironment.  why?
    class ParallelMachineFinalize
    {
      bool m_need_to_finalize;
    public:
      ParallelMachineFinalize(bool need_to_finalize=false) : m_need_to_finalize(need_to_finalize) {}
      ~ParallelMachineFinalize()
      {
        if (m_need_to_finalize)
          {
            stk::parallel_machine_finalize();
          }
      }
    };

    class RunEnvironment : public ParallelMachineFinalize
    {

    public:
      // Will initialize a comm
      RunEnvironment(int *argc, char ***argv, bool debug=false);

      // Assumes already-initialized comm
      RunEnvironment(int *argc, char ***argv, stk::ParallelMachine comm, bool debug=false);

      //int processCommandLine() { return processCommandLine(m_argc, m_argv); }

      ~RunEnvironment();

      std::string
      build_log_description(const std::string &           working_directory,
                            int                           parallel_rank,
                            int                           parallel_size);

      int get_argc() { return m_argc; }
      char **get_argv() { return m_argv; }



      // command line options
      Teuchos::CommandLineProcessor clp;

      std::string output_log_opt;
      std::string dw_opt;
      std::string timer_opt;
      std::string directory_opt;

      std::string pout_opt;
      std::string dout_opt;
      std::string runtest_opt;

      int help_opt;

      // data
      const stk::ParallelMachine    m_comm;
      static std::string            m_workingDirectory;

    private:

      bool                          m_need_to_finalize;
      bool                          m_debug;
      bool                          m_processCommandLine_invoked;
      std::string                  *m_argv_new;
      int                           m_argc;
      char                        **m_argv;

      //ParallelMachineFinalize       m_par_finalize;

      // shared constructor implementation; do not call directly
      void internal_initialize(int argc, char **argv);
      void bootstrap();

      void setSierraOpts(int procRank, int argc, char* argv[]);
    };

    int processCommandLine(Teuchos::CommandLineProcessor& clp_in, int argc, char **argv);

    int processCLP(Teuchos::CommandLineProcessor& clp_in, int procRank, int argc, char* argv[]);

    std::string get_working_directory();

    int setFileNames(std::string& fullmesh, std::string& meshFileName, std::string& errString);

    void runCommand(std::string command);

    void doLoadBalance(stk::ParallelMachine comm, std::string meshFileName);

    void printHelp(Teuchos::CommandLineProcessor& clp_in);

  } // namespace percept

#endif // percept_RunEnvironment_hpp