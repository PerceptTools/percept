#ifndef test_structs_hpp
#define test_structs_hpp

#include <Kokkos_Core.hpp>
#include <percept/MeshType.hpp>
#include <array>

#include <percept/MeshType.hpp>
#include <percept/structured/StructuredBlock.hpp>
#include <Kokkos_Vector.hpp>

#include <percept/PerceptMesh.hpp>
#include <percept/structured/BlockStructuredGrid.hpp>

/////////////////////////////////////////////////////////////////////////////////////////////
      void failsKokkos()
      {//doesn't even get called, which is good I suppose
          std::cout << "You can't do this from GPU\n";
      }

      KOKKOS_INLINE_FUNCTION
      double det_test(double m[3][3])
      {
          printf("det_test call\n");
          return m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
                    + m[0][1]*(m[2][0]*m[1][2] - m[1][0]*m[2][2])
                    + m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]);
      }

      KOKKOS_INLINE_FUNCTION
      bool jacobian_matrix_3D_test(double &detJ, const  double * x0, const double * x1, const double * x2, const double * x3)
      {
          printf("jacobian_matrix_3D_test call\n");
          double  A[3][3];

          A[0][0] = (x1[0] - x0[0]);
          A[0][1] = (x2[0] - x0[0]);
          A[0][2] = (x3[0] - x0[0]);

          A[1][0] = (x1[1] - x0[1]);
          A[1][1] = (x2[1] - x0[1]);
          A[1][2] = (x3[1] - x0[1]);

          A[2][0] = (x1[2] - x0[2]);
          A[2][1] = (x2[2] - x0[2]);
          A[2][2] = (x3[2] - x0[2]);

          detJ = det_test(A);
          return detJ < 0.0;
      }

      KOKKOS_INLINE_FUNCTION
      bool SGridJacobianUtil_test(double& averageJ, double detJ[8])
      {
          printf("calling SGrisjacobianUtil_test\n");
        const int locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                4, 2 }, { 7, 5, 6, 3 } };

        bool metric_invalid = false;

        //const int A0 = sgrid->m_access_ordering[0], A1 = sgrid->m_access_ordering[1], A2 = sgrid->m_access_ordering[2];
//        const int A0 = 0, A1 = 1, A2 = 2;

        double v_i[8][3];
//        unsigned indx[3] = { 0, 0, 0 };
//        unsigned II[3] = { 0, 0, 0 };

//        unsigned cnt = 0;
//          for (indx[2] = 0; indx[2] < 2; ++indx[2]) {
//              II[2] = indx[2] + cell_ijkb[2];
//              for (indx[1] = 0; indx[1] < 2; ++indx[1]) {
//                  II[1] = indx[1] + cell_ijkb[1];
//                  for (indx[0] = 0; indx[0] < 2; ++indx[0]) {
//                      II[0] = indx[0] + cell_ijkb[0];
//                      for (unsigned ic = 0; ic < 3; ++ic) {
//                          v_i[cnt][ic] = coords(II[A0], II[A1], II[A2], ic);
//                      }
//                      ++cnt;
//                  }
//              }
//          }

        for (int i = 0; i < 8; ++i) {
            bool mi = jacobian_matrix_3D_test(
                    detJ[i],
                    v_i[locs_hex_dev[i][0]], v_i[locs_hex_dev[i][1]],
                    v_i[locs_hex_dev[i][2]], v_i[locs_hex_dev[i][3]]);

            metric_invalid = metric_invalid || mi;
        }

        {
            averageJ = 0.0;
            for (unsigned i = 0; i < 8; i++)
                averageJ += detJ[i];
            averageJ /= double(8);
        }

        return metric_invalid;
      }

      struct testMemberStruct
      {
          const double m_beta_mult = 0.5;
          KOKKOS_INLINE_FUNCTION
          double metric_test(bool& valid) const
          {
              printf("calling metric_test\n");
               valid = true;

                double A_ = 0.0, W_ = 0.0; // current and reference detJ
                double nodal_A[8], nodal_W[8];

                SGridJacobianUtil_test(A_, nodal_A);
                SGridJacobianUtil_test(W_, nodal_W);
                double val_untangle=0.0;

                for (int i=0; i < 8; i++)
                  {
                    double detAi = nodal_A[i];
                    double detWi = nodal_W[i];

                    if (detAi <= 0.) valid = false;

                    double vv = m_beta_mult*detWi - detAi;

                    (vv<0.0 ? vv=0.0 : vv);
                    val_untangle += vv*vv;
                  }
                return val_untangle;
          }
      };

      struct testStruct
      {
          mutable testMemberStruct tms;

          void run()
          {
              double total;
              Kokkos::parallel_reduce(20,*this,total);
              std::cout << "The total is ... "<<total << std::endl;
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index, double& mtot_loc) const
          {
            const_cast<testStruct *>(this)->operator()(index, mtot_loc);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index, double& mtot_loc)
          {
              printf("Operator call\n");
              bool somebool;
              mtot_loc += tms.metric_test(somebool);
//            failsKokkos();//somehow it just ignores this call
          }
      };

      struct datatypetest_doubleParallelReduce
      {
        double commonval;
        double total;
        unsigned sz;
        Kokkos::View<double*, percept::DataLayout, percept::MemSpace> myview;
        Kokkos::View<double*, percept::DataLayout, percept::MemSpace>::HostMirror mymirror;

        datatypetest_doubleParallelReduce(unsigned sz_in, double divider)
        {
            if(sz_in==0)
                sz_in=1;

            total = 0.0;
            commonval = 1.0/divider;
            sz = sz_in;
            Kokkos::View<double*, percept::DataLayout, percept::MemSpace> inter("myview",sz);
            myview=inter;

            mymirror = Kokkos::create_mirror_view(myview);
            for(unsigned i=0;i<sz;i++)
                mymirror[i]=commonval;
            Kokkos::deep_copy(myview,mymirror);
        }

        void reduce()
        {
            double local = total;
            Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0,sz),*this,local);
            total = local;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned& index, double& local_sum) const
        {
          const_cast<datatypetest_doubleParallelReduce *>(this)->operator()(index, local_sum);
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned& index, double& local_sum)
        {
            local_sum += myview(index);
        }
      };

      struct datatypetest_long_doubleParallelReduce
      {//it appears long doubles cause illegal memory access errors when summed into on a GPU using CUDA. This is due to the fact that CUDA GPU's only use single or double floating pointer precision
        double commonval;
        long double total;
        unsigned sz;
        Kokkos::View<double*, percept::DataLayout, percept::MemSpace> myview;
        Kokkos::View<double*, percept::DataLayout, percept::MemSpace>::HostMirror mymirror;

        datatypetest_long_doubleParallelReduce(unsigned sz_in, double divider)
        {
            if(sz_in==0)
                sz_in=1;

            total = 0.0;
            commonval = 1.0/divider;
            sz = sz_in;
            Kokkos::View<double*, percept::DataLayout, percept::MemSpace> inter("myview",sz);
            myview=inter;

            mymirror = Kokkos::create_mirror_view(myview);
            for(unsigned i=0;i<sz;i++)
                mymirror[i]=commonval;
            Kokkos::deep_copy(myview,mymirror);
        }

        void reduce()
        {
            long double local = total;
            Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0,sz),*this,local);
            total = local;
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned& index, long double& local_sum) const
        {
          const_cast<datatypetest_long_doubleParallelReduce *>(this)->operator()(index, local_sum);
        }

        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned& index, long double& local_sum)
        {
            local_sum += myview(index);
        }
      };


      class virtual_functions_kokkos_test
      {
      public:
          KOKKOS_INLINE_FUNCTION
          virtual void doSomething(unsigned arg) const = 0;

//        virtual ~virtual_functions_kokkos_test() = 0;
      };

      class inherit_virtual_test : public virtual_functions_kokkos_test
      {
      public:
          //It appears you can build a class inheriting from a pure virtual class using Kokkos::Cuda. However, the override on the virtual function doesn't appear to callable from the operator (it can be called from run function).
          //The code will compile and run with the virtual function called in the operator, but its output never displays from the operator

          KOKKOS_INLINE_FUNCTION
          virtual void doSomething(unsigned arg) const {printf("a thing was done virt ... %d\n",arg);}

          KOKKOS_INLINE_FUNCTION
          void doSomething_non_virt(unsigned arg) const {printf("a thing was done non_virt ... %d\n",arg);}


          void run()
          {
              doSomething(99);
              Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,50),*this);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index) const
          {
            const_cast<inherit_virtual_test *>(this)->operator()(index);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index)
          {
              printf("Operator\n");
              doSomething_non_virt(index); //this appears to get called. However, it doesn't seem to get called very time; not all indices are printed
              doSomething(index); //doesn't even seem to get called
          }
      };

      class no_virtual_test
      {//this class functions as expected
      public:
          KOKKOS_INLINE_FUNCTION
          void doSomething_non_virt(unsigned arg) const {printf("a thing was done non_virt nvt... %d\n",arg);}


          void run()
          {
              Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,50),*this);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index) const
          {
            const_cast<no_virtual_test *>(this)->operator()(index);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index)
          {
              printf("Operator nvt\n");
              doSomething_non_virt(index);
          }
      };


      class using_std_ds_cuda
      {//If built and run with cuda, this triggers an illegal memory access error. This is curious, as it begs the question why did this never happen with the structured refinement code?
          //What is the essential difference between the refinement and this from a memory standpoint?
          Kokkos::View<int*, percept::DataLayout,percept::MemSpace> myview;
          Kokkos::View<int*, percept::DataLayout,percept::MemSpace>::HostMirror mymirror;
      public:
          using_std_ds_cuda()
          {
              Kokkos::View<int*, percept::DataLayout,percept::MemSpace> interim("myview",10);
              myview = interim;
              mymirror = Kokkos::create_mirror_view(myview);
              for(int i=0;i<10;i++)
                  mymirror(i)=0;
              Kokkos::deep_copy(myview,mymirror);
          }

          void run()
          {
              printf("About to run PF  using_std_ds_cuda  \n");
              Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,10),*this);
              Kokkos::deep_copy(mymirror,myview);
              bool valid = true;
              for(int i=0;i<10;i++){
                  if(mymirror(i)!=5){
                      printf("incorrect value at index %d\n",i);
                      valid = false;
                      break;
                  }
              }
              if(valid)
                  printf("executed as intended using_std_ds_cuda\n\n\n");
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index) const
          {
            const_cast<using_std_ds_cuda *>(this)->operator()(index);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index)
          {
              printf("Operator call at index %d   using_std_ds_cuda\n",index);
              std::array<int,1> arr;
              arr[0]=5;
              printf("Set the stda::array value at index %d\n",index);//never reaches here
              myview(index)=arr[0];
              printf("set the view's value at index %d\n",index);
          }
      };




      class using_kokkos_ds_cuda
      {
          Kokkos::View<int*, percept::DataLayout,percept::MemSpace> myview;
          Kokkos::View<int*, percept::DataLayout,percept::MemSpace>::HostMirror mymirror;
      public:
          using_kokkos_ds_cuda()
          {
              Kokkos::View<int*, percept::DataLayout,percept::MemSpace> interim("myview",10);
              myview = interim;
              mymirror = Kokkos::create_mirror_view(myview);
              for(int i=0;i<10;i++)
                  mymirror(i)=0;
              Kokkos::deep_copy(myview,mymirror);
          }

          void run()
          {
              printf("About to run PF  using_kokkos_ds_cuda  \n");
              Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,10),*this);
              Kokkos::deep_copy(mymirror,myview);
              bool valid = true;
              for(int i=0;i<10;i++){
                  if(mymirror(i)!=5){
                      printf("incorrect value at index %d\n",i);
                      valid = false;
                      break;
                  }
              }
              if(valid)
                  printf("executed as intended using_kokkos_ds_cuda\n\n\n");
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index) const
          {
            const_cast<using_kokkos_ds_cuda *>(this)->operator()(index);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const unsigned& index)
          {
              printf("Operator call at index %d   using_kokkos_ds_cuda\n",index);
              Kokkos::Array<int,1> arr;
              arr[0]=5;
              printf("Set the stda::array value at index %d\n",index);//never reaches here
              myview(index)=arr[0];
              printf("set the view's value at index %d\n",index);
          }
      };

      class release_vs_debug_openmp_fail_reduce {
private:
    const unsigned indxval = 6;
    const unsigned div = 2;
    Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace> myview;
public:
    release_vs_debug_openmp_fail_reduce(unsigned p_size) {
        Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace> vals(
                "myview", p_size);
        myview = vals;
        Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace>::HostMirror mir;
        mir = Kokkos::create_mirror_view(myview);
        for (unsigned i = 0; i < p_size; i++) {
            mir(i) = indxval;
        }
        Kokkos::deep_copy(myview, mir);
    }

    unsigned reduce() {
        unsigned total;
        Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0,myview.size()),*this,total);
        return total;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index, unsigned& mtot_loc) const {
        const_cast<release_vs_debug_openmp_fail_reduce *>(this)->operator()(index,mtot_loc);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index, unsigned& mtot_loc) {
        mtot_loc += (myview(index)/div);
    }

};

      class release_vs_debug_openmp_fail_for {
private:
    const unsigned indxval = 6;
    const unsigned div = 2;
    Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace>::HostMirror mir;
    Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace> myview;
public:
    release_vs_debug_openmp_fail_for(unsigned p_size) {
        Kokkos::View<unsigned*, percept::DataLayout, percept::MemSpace> vals(
                "myview", p_size);
        myview = vals;
        mir = Kokkos::create_mirror_view(myview);
    }

    bool fill(){
        Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,myview.size()),*this);
        Kokkos::deep_copy(mir,myview);
        for(unsigned i=0;i<mir.size();i++)
        {
            if (mir(i)!=indxval/div)
            {//if a value wasn't set correctly
                return false;
            }
        }
        return true;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) const {
        const_cast<release_vs_debug_openmp_fail_for *>(this)->operator()(index);

    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned& index) {
        myview(index)=(indxval/div);
    }

};



      struct det_arr_tester
      {
          uint64_t  iteration_sz;
          det_arr_tester(uint64_t its_in)
          {
              iteration_sz=its_in;
          }

          void run()
          {
              double tot=0;
              Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0,iteration_sz),*this,tot);


          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const uint64_t & index, double& loc_tot) const
          {
            const_cast<det_arr_tester *>(this)->operator()(index, loc_tot);
          }

          KOKKOS_INLINE_FUNCTION
          bool set(double & det, double m[3][3]) const
          {
              m[0][0]=.5;
              m[0][1]=2.5;
              m[0][2]=5.0;

              m[1][0]=3.3;
              m[1][1]=7.9;
              m[1][2]=12.2;

              m[2][0]=47.8;
              m[2][1]=-100.8;
              m[2][2]=-1.6;

              det = Det(m);

              return det < 0.0;
          }

          KOKKOS_INLINE_FUNCTION
          double Det(double m[3][3]) const
          {
              return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                     + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                     + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
          }

          KOKKOS_INLINE_FUNCTION
          void operator()(const uint64_t & index, double& loc_tot)
          {

              double det;
              double m[3][3];

              set(det,m);

              loc_tot +=det;

          }
      };

  /*struct*/class det_arr_tester_w_data {
  private:
      Kokkos::View<double**, percept::DataLayout, percept::MemSpace> mView;

      uint64_t m_sz;
  public:
      det_arr_tester_w_data(uint64_t sz) {
          Kokkos::View<double**, percept::DataLayout, percept::MemSpace> viewin("mView", 8, sz);
          mView = viewin;

          Kokkos::View<double**, percept::DataLayout, percept::MemSpace>::HostMirror mir;
          mir = Kokkos::create_mirror_view(mView);

          for (unsigned iElem = 0; iElem < sz; iElem++) {
              mir(0, iElem) = (double) iElem + 1.5;
              mir(1, iElem) = (double) iElem * .5 + 2.9;
              mir(2, iElem) = (double) iElem * .24 + 3.7;
              mir(3, iElem) = (double) iElem * 1.5 + 4.14545;
              mir(4, iElem) = (double) iElem * 2.6 + 5.15;
              mir(5, iElem) = (double) iElem * .0005 + 6.22;
              mir(6, iElem) = (double) iElem * .67 + 7.3;
              mir(7, iElem) = (double) iElem * (double) (-1.0 / iElem) + 8.4;
          }
          Kokkos::deep_copy(mView, mir);
          m_sz = sz;
      }

      void run() {
          double tot = 0;
          Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0, m_sz), *this,
                  tot);
      }


      KOKKOS_INLINE_FUNCTION
      void operator()(const uint64_t & index, double& loc_tot) const {
          const_cast<det_arr_tester_w_data *>(this)->operator()(index, loc_tot);
      }
  private:
      KOKKOS_INLINE_FUNCTION
      bool set(double & det, double m[3][3], double coords[8]) const {
          m[0][0] = coords[0] + coords[1];
          m[0][1] = coords[1] + coords[2];
          m[0][2] = coords[2] + coords[3];

          m[1][0] = coords[3] + coords[4];
          m[1][1] = coords[4] + coords[5];
          m[1][2] = coords[5] + coords[6];

          m[2][0] = coords[6] + coords[7];
          m[2][1] = coords[7] + coords[0];
          m[2][2] = coords[7] + coords[1];

          det = Det(m);

          return det < 0.0;
      }

      KOKKOS_INLINE_FUNCTION
      double Det(double m[3][3]) const {
          return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                  + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                  + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const uint64_t & index, double& loc_tot) {

          double coords[8];
          double m[3][3];
          double det;
          for (unsigned iCoord = 0; iCoord < 8; iCoord++)
              coords[iCoord] = mView(iCoord, index);
          set(det, m, coords);

          loc_tot += det;

      }

  };

  class m_functions
  {
  private:
      Kokkos::View<double**, percept::DataLayout, percept::MemSpace> mView;
  public:
      m_functions(Kokkos::View<double**, percept::DataLayout, percept::MemSpace> viewin) { mView = viewin; }

      KOKKOS_INLINE_FUNCTION
      double fake_metric(double coords[8])
      {
          double t;
          SGJU(t,coords);
          SGJU(t,coords);
          return t;
      }

  private:

      KOKKOS_INLINE_FUNCTION
      bool SGJU(double& in, double coords[8])
      {
          double m[3][3];
          for(int i=0;i<8;i++)
              set(in,m,coords);
          return in < 0.0;
      }

      KOKKOS_INLINE_FUNCTION
      bool set(double & det, double m[3][3], double coords[8]) const {
          m[0][0] = coords[0] + coords[1];
          m[0][1] = coords[1] + coords[2];
          m[0][2] = coords[2] + coords[3];

          m[1][0] = coords[3] + coords[4];
          m[1][1] = coords[4] + coords[5];
          m[1][2] = coords[5] + coords[6];

          m[2][0] = coords[6] + coords[7];
          m[2][1] = coords[7] + coords[0];
          m[2][2] = coords[7] + coords[1];

          det = Det(m);

          return det < 0.0;
      }


      KOKKOS_INLINE_FUNCTION
      double Det(double m[3][3]) const {
          return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                  + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                  + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
      }
  };

  struct function_user {

      Kokkos::View<double**, percept::DataLayout, percept::MemSpace> mView;
      m_functions mfunct;
      uint64_t m_sz;

      function_user(uint64_t sz) : mfunct(mView) {
          Kokkos::View<double**, percept::DataLayout, percept::MemSpace> viewin("mView", 8, sz);
          mView = viewin;

          Kokkos::View<double**, percept::DataLayout, percept::MemSpace>::HostMirror mir;
          mir = Kokkos::create_mirror_view(mView);

          for (unsigned iElem = 0; iElem < sz; iElem++) {
              mir(0, iElem) = (double) iElem + 1.5;
              mir(1, iElem) = (double) iElem * .5 + 2.9;
              mir(2, iElem) = (double) iElem * .24 + 3.7;
              mir(3, iElem) = (double) iElem * 1.5 + 4.14545;
              mir(4, iElem) = (double) iElem * 2.6 + 5.15;
              mir(5, iElem) = (double) iElem * .0005 + 6.22;
              mir(6, iElem) = (double) iElem * .67 + 7.3;
              mir(7, iElem) = (double) iElem * (double) (-1.0 / iElem) + 8.4;
          }
          Kokkos::deep_copy(mView, mir);
          m_sz = sz;

      }

      void run() {
          double tot = 0;
          Kokkos::parallel_reduce(Kokkos::RangePolicy<percept::ExecSpace>(0, m_sz), *this,
                  tot);
      }


      KOKKOS_INLINE_FUNCTION
      void operator()(const uint64_t & index, double& loc_tot) const {
          const_cast<function_user *>(this)->operator()(index, loc_tot);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const uint64_t & index, double& loc_tot) {

          double coords[8];
  //        double m[3][3];
          double det;
          for (unsigned iCoord = 0; iCoord < 8; iCoord++)
              coords[iCoord] = mView(iCoord, index);
          det = mfunct.fake_metric(coords);

          loc_tot += det;
      }

  };


  struct test_gpu_scalings_mesh {

    double m_beta_mult;

    bool m_use_ref_mesh;
    bool m_untangling;

    percept::StructuredGrid::MTField *m_coord_field_current;
    percept::StructuredGrid::MTField *m_coord_field_original;

    percept::StructuredGrid::MTField::Array4D m_coords_current_iterator;
    percept::StructuredGrid::MTField::Array4D m_coords_original_iterator;
    Kokkos::View<unsigned**, percept::DataLayout, percept::MemSpace> logical_elems;
    unsigned nele;

    test_gpu_scalings_mesh(percept::PerceptMesh *eMesh) :
            m_beta_mult(0.05) {
        if (eMesh->get_block_structured_grid()) {
            std::shared_ptr<percept::BlockStructuredGrid> bsg =
                    eMesh->get_block_structured_grid();

            m_coord_field_current = bsg->m_fields["coordinates"].get();
            m_coord_field_original = bsg->m_fields["coordinates_NM1"].get();
            m_coords_current_iterator =
                    *m_coord_field_current->m_block_fields[0]; //orient iterators at the begining of blocks
            m_coords_original_iterator =
                    *m_coord_field_original->m_block_fields[0];
            m_use_ref_mesh = true; //will get changed by the smoother
            m_untangling = true; //will get changed by the smoother

            std::vector<percept::StructuredCellIndex> elems_from_block;
            bsg->get_elements_of_sb(elems_from_block, 0);

            Kokkos::View<unsigned**, percept::DataLayout, percept::MemSpace>::HostMirror elems_mirror;

            //populate the element ijk's of a particular block
            Kokkos::resize(logical_elems, elems_from_block.size(), 3);
            elems_mirror = Kokkos::create_mirror_view(logical_elems);
            bsg->get_elements_of_sb(elems_from_block, 0);
            for (unsigned iElem = 0; iElem < elems_from_block.size(); iElem++) {
                elems_mirror(iElem, 0) = elems_from_block[iElem][0];
                elems_mirror(iElem, 1) = elems_from_block[iElem][1];
                elems_mirror(iElem, 2) = elems_from_block[iElem][2];
            }
            Kokkos::deep_copy(logical_elems, elems_mirror); //sync back to device

            nele = elems_from_block.size();
        }
    }

    KOKKOS_INLINE_FUNCTION
    double det(double m[3][3]) const {
        //      printf("calling member det\n");

        return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
    } //double det

    KOKKOS_INLINE_FUNCTION
    bool jacobian_matrix_3D(double &detJ, double A[3][3], const double * x0,
            const double * x1, const double * x2, const double * x3) const {
        //      printf("calling member jacobian_matrix_3D\n");
        A[0][0] = (x1[0] - x0[0]);
        A[0][1] = (x2[0] - x0[0]);
        A[0][2] = (x3[0] - x0[0]);

        A[1][0] = (x1[1] - x0[1]);
        A[1][1] = (x2[1] - x0[1]);
        A[1][2] = (x3[1] - x0[1]);

        A[2][0] = (x1[2] - x0[2]);
        A[2][1] = (x2[2] - x0[2]);
        A[2][2] = (x3[2] - x0[2]);

        detJ = det(A);

        return detJ < 0.0;
    } //double jacobian_matrix_3D

    KOKKOS_INLINE_FUNCTION
    bool sGridJacobianUtil(double detJ[8],
    /*percept::StructuredGrid::MTField::Array4D coords,
     Kokkos::Array<unsigned, 3> cell_ijk*/double coords[8][3],
            double J[8][3][3]) const {

        const int locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                4, 2 }, { 7, 5, 6, 3 } };

        bool metric_invalid = false;

        for (int i = 0; i < 8; ++i) {
            bool mi = jacobian_matrix_3D(detJ[i], J[i],
                    coords[locs_hex_dev[i][0]], coords[locs_hex_dev[i][1]],
                    coords[locs_hex_dev[i][2]], coords[locs_hex_dev[i][3]]);

            metric_invalid = metric_invalid || mi;
        }

        return metric_invalid;
    } //bool sGriJacobianUtil

    KOKKOS_INLINE_FUNCTION
    double metric(/*Kokkos::Array<unsigned, 3> element*/double elem_curr[8][3],
            double elem_org[8][3], bool& valid) const {

        valid = true;
        double nodal_A[8], nodal_W[8];
        double J_A[8][3][3];
        double J_W[8][3][3];

        sGridJacobianUtil(/*A_,*/nodal_A,/* m_coords_current_iterator, element*/
                elem_curr, J_A);
        sGridJacobianUtil(/*W_,*/nodal_W,/* m_coords_original_iterator, element*/
                elem_org, J_W);

        if (m_untangling) {
            double val_untangle = 0.0;
            for (int i = 0; i < 8; i++) {

                double detAi = nodal_A[i];
                double detWi = nodal_W[i];

                if (detAi <= 0.)
                    valid = false;

                double vv = m_beta_mult * detWi - detAi;
                (vv < 0.0 ? vv = 0.0 : vv);
                val_untangle += vv * vv;
            }
            return val_untangle;
        }

        return 0.0;
    } //metric

    void run() {
        double tot = 0;
        Kokkos::parallel_reduce(
                Kokkos::RangePolicy<percept::ExecSpace>(0, nele), *this, tot);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned & index, double& loc_tot) const {
        const_cast<test_gpu_scalings_mesh *>(this)->operator()(index, loc_tot);
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(const unsigned & index, double& loc_tot) {
        bool local_valid = false;

        double v_i_current[8][3];
        double v_i_org[8][3];
        unsigned indx[3] = { 0, 0, 0 };
        unsigned II[3] = { 0, 0, 0 };
        unsigned cell_ijk[3] = { logical_elems(index, 0), logical_elems(index,
                1), logical_elems(index, 2) };
        const int A0 = 0, A1 = 1, A2 = 2;

        unsigned cnt = 0;
        for (indx[2] = 0; indx[2] < 2; ++indx[2]) {
            II[2] = indx[2] + cell_ijk[2];
            for (indx[1] = 0; indx[1] < 2; ++indx[1]) {
                II[1] = indx[1] + cell_ijk[1];
                for (indx[0] = 0; indx[0] < 2; ++indx[0]) {
                    II[0] = indx[0] + cell_ijk[0];
                    for (unsigned ic = 0; ic < 3; ++ic) {
                        v_i_current[cnt][ic] = m_coords_current_iterator(II[A0],
                                II[A1], II[A2], ic);
                        v_i_org[cnt][ic] = m_coords_original_iterator(II[A0],
                                II[A1], II[A2], ic);
                    }
                    ++cnt;
                }
            }
        }

        percept::Double mm = metric(v_i_current, v_i_org, local_valid);

        loc_tot += mm;
    }
};


  struct test_gpu_scalings_functions_only
  {
      unsigned nele;

      test_gpu_scalings_functions_only(unsigned nelein) : nele(nelein) {}

      KOKKOS_INLINE_FUNCTION
      double det(double m[3][3]) const {
          //      printf("calling member det\n");

          return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                  + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                  + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
      } //double det

      KOKKOS_INLINE_FUNCTION
      bool jacobian_matrix_3D(double &detJ, double A[3][3], const double * x0,
              const double * x1, const double * x2, const double * x3) const {
          //      printf("calling member jacobian_matrix_3D\n");
          A[0][0] = (x1[0] - x0[0]);
          A[0][1] = (x2[0] - x0[0]);
          A[0][2] = (x3[0] - x0[0]);

          A[1][0] = (x1[1] - x0[1]);
          A[1][1] = (x2[1] - x0[1]);
          A[1][2] = (x3[1] - x0[1]);

          A[2][0] = (x1[2] - x0[2]);
          A[2][1] = (x2[2] - x0[2]);
          A[2][2] = (x3[2] - x0[2]);

          detJ = det(A);

          return detJ < 0.0;
      } //double jacobian_matrix_3D

      KOKKOS_INLINE_FUNCTION
      bool sGridJacobianUtil(double detJ[8],
      /*percept::StructuredGrid::MTField::Array4D coords,
       Kokkos::Array<unsigned, 3> cell_ijk*/double coords[8][3],
              double J[8][3][3]) const {

          const int locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                  3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                  4, 2 }, { 7, 5, 6, 3 } };

          bool metric_invalid = false;

          for (int i = 0; i < 8; ++i) {
              bool mi = jacobian_matrix_3D(detJ[i], J[i],
                      coords[locs_hex_dev[i][0]], coords[locs_hex_dev[i][1]],
                      coords[locs_hex_dev[i][2]], coords[locs_hex_dev[i][3]]);

              metric_invalid = metric_invalid || mi;
          }

          return metric_invalid;
      } //bool sGriJacobianUtil

      KOKKOS_INLINE_FUNCTION
      double metric(/*Kokkos::Array<unsigned, 3> element*/double elem_curr[8][3],
              double elem_org[8][3], bool& valid) const {

          valid = true;
          double nodal_A[8], nodal_W[8];
          double J_A[8][3][3];
          double J_W[8][3][3];

          sGridJacobianUtil(/*A_,*/nodal_A,/* m_coords_current_iterator, element*/
                  elem_curr, J_A);
          sGridJacobianUtil(/*W_,*/nodal_W,/* m_coords_original_iterator, element*/
                  elem_org, J_W);

          if (true) {
              double val_untangle = 0.0;
              for (int i = 0; i < 8; i++) {

                  double detAi = nodal_A[i];
                  double detWi = nodal_W[i];

                  if (detAi <= 0.)
                      valid = false;

                  double vv = 0.05 * detWi - detAi;
                  (vv < 0.0 ? vv = 0.0 : vv);
                  val_untangle += vv * vv;
              }
              return val_untangle;
          }

          return 0.0;
      } //metric

      void run() {
          double tot = 0;
          Kokkos::parallel_reduce(
                  Kokkos::RangePolicy<percept::ExecSpace>(0, nele), *this, tot);
          std::cout << "The total is .. " << tot << std::endl;
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned & index, double& loc_tot) const {
          const_cast<test_gpu_scalings_functions_only *>(this)->operator()(index, loc_tot);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const unsigned & index, double& loc_tot) {
          bool local_valid = false;

          double v_i_current[8][3];
          double v_i_org[8][3];


          for(unsigned i=0;i<8;i++)
              for(unsigned j=0;j<3;j++)
              {
                  v_i_current[i][j] = 1.5*(i+j);
                  v_i_org[i][j] = 1.5*(i-j);
              }

          percept::Double mm = metric(v_i_current, v_i_org, local_valid);

          loc_tot += mm;
      }
  };

  struct test_gpu_scalings_custom_data
   {
       unsigned nele;
       Kokkos::View<double**, percept::DataLayout,percept::MemSpace> coordsView;

       test_gpu_scalings_custom_data(unsigned nelein) : nele(nelein)
       {
           Kokkos::View<double**, percept::DataLayout,percept::MemSpace> vin("coordsView",nele,3);
           coordsView=vin;
           Kokkos::View<double**, percept::DataLayout,percept::MemSpace>::HostMirror mir = Kokkos::create_mirror_view(coordsView);
           for(unsigned iElem=0;iElem<nele;iElem++)
                   for(unsigned iCoord=0;iCoord<3;iCoord++)
                       mir(iElem,iCoord) = (iCoord + .5)*(iElem+.0001);
           Kokkos::deep_copy(coordsView,mir);
       }

       KOKKOS_INLINE_FUNCTION
       double det(double m[3][3]) const {
           //      printf("calling member det\n");

           return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
                   + m[0][1] * (m[2][0] * m[1][2] - m[1][0] * m[2][2])
                   + m[0][2] * (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
       } //double det

       KOKKOS_INLINE_FUNCTION
       bool jacobian_matrix_3D(double &detJ, double A[3][3], const double * x0,
               const double * x1, const double * x2, const double * x3) const {
           //      printf("calling member jacobian_matrix_3D\n");
           A[0][0] = (x1[0] - x0[0]);
           A[0][1] = (x2[0] - x0[0]);
           A[0][2] = (x3[0] - x0[0]);

           A[1][0] = (x1[1] - x0[1]);
           A[1][1] = (x2[1] - x0[1]);
           A[1][2] = (x3[1] - x0[1]);

           A[2][0] = (x1[2] - x0[2]);
           A[2][1] = (x2[2] - x0[2]);
           A[2][2] = (x3[2] - x0[2]);

           detJ = det(A);

           return detJ < 0.0;
       } //double jacobian_matrix_3D

       KOKKOS_INLINE_FUNCTION
       bool sGridJacobianUtil(double detJ[8],
       /*percept::StructuredGrid::MTField::Array4D coords,
        Kokkos::Array<unsigned, 3> cell_ijk*/double coords[8][3],
               double J[8][3][3]) const {

           const int locs_hex_dev[8][4] = { { 0, 1, 2, 4 }, { 1, 3, 0, 5 }, { 2, 0,
                   3, 6 }, { 3, 2, 1, 7 }, { 4, 6, 5, 0 }, { 5, 4, 7, 1 }, { 6, 7,
                   4, 2 }, { 7, 5, 6, 3 } };

           bool metric_invalid = false;

           for (int i = 0; i < 8; ++i) {
               bool mi = jacobian_matrix_3D(detJ[i], J[i],
                       coords[locs_hex_dev[i][0]], coords[locs_hex_dev[i][1]],
                       coords[locs_hex_dev[i][2]], coords[locs_hex_dev[i][3]]);

               metric_invalid = metric_invalid || mi;
           }

           return metric_invalid;
       } //bool sGriJacobianUtil

       KOKKOS_INLINE_FUNCTION
       double metric(/*Kokkos::Array<unsigned, 3> element*/double elem_curr[8][3],
               double elem_org[8][3], bool& valid) const {

           valid = true;
           double nodal_A[8], nodal_W[8];
           double J_A[8][3][3];
           double J_W[8][3][3];

           sGridJacobianUtil(/*A_,*/nodal_A,/* m_coords_current_iterator, element*/
                   elem_curr, J_A);
           sGridJacobianUtil(/*W_,*/nodal_W,/* m_coords_original_iterator, element*/
                   elem_org, J_W);

           if (true) {
               double val_untangle = 0.0;
               for (int i = 0; i < 8; i++) {

                   double detAi = nodal_A[i];
                   double detWi = nodal_W[i];

                   if (detAi <= 0.)
                       valid = false;

                   double vv = 0.05 * detWi - detAi;
                   (vv < 0.0 ? vv = 0.0 : vv);
                   val_untangle += vv * vv;
               }
               return val_untangle;
           }

           return 0.0;
       } //metric

       void run() {
           double tot = 0;
           Kokkos::parallel_reduce(
                   Kokkos::RangePolicy<percept::ExecSpace>(0, nele), *this, tot);
           std::cout << "The total is .. " << tot << std::endl;
       }

       KOKKOS_INLINE_FUNCTION
       void operator()(const unsigned & index, double& loc_tot) const {
           const_cast<test_gpu_scalings_custom_data *>(this)->operator()(index, loc_tot);
       }

       KOKKOS_INLINE_FUNCTION
       void operator()(const unsigned & index, double& loc_tot) {
           bool local_valid = false;

           double v_i_current[8][3];
           double v_i_org[8][3];


           for(unsigned iNode=0;iNode<8;iNode++)
               for(unsigned jCoord=0;jCoord<3;jCoord++)
               {
                   v_i_current[iNode][jCoord] = (iNode+1.0)*4*coordsView(index,jCoord);
                   v_i_org[iNode][jCoord] = (iNode+1.0)*coordsView(index,jCoord);
               }

           percept::Double mm = metric(v_i_current, v_i_org, local_valid);

           loc_tot += mm;
       }
   };

  struct lambda_struct
  {
      Kokkos::View<double*,percept::DataLayout,percept::MemSpace> v;

      lambda_struct()
      {
          Kokkos::View<double*,percept::DataLayout,percept::MemSpace> vin("v",100);
          v=vin;
      }

      void set() {

        auto sv = subview(v,Kokkos::make_pair(0,50));

        Kokkos::parallel_for(Kokkos::RangePolicy<percept::ExecSpace>(0,sv.size()), KOKKOS_LAMBDA (unsigned index) {
            sv[index] = 7.0;
    });
    }
  };


#endif
