// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIexpdIdunedIappdIusersdIcnguyendIlblpwgtoolsdICAFAnadItutedIcpv_joint_Multi_noSystem_C_ACLiC_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/exp/dune/app/users/cnguyen/lblpwgtools/CAFAna/tute/./cpv_joint_Multi_noSystem.C"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *VectorCompareResult_Dictionary();
   static void VectorCompareResult_TClassManip(TClass*);
   static void *new_VectorCompareResult(void *p = 0);
   static void *newArray_VectorCompareResult(Long_t size, void *p);
   static void delete_VectorCompareResult(void *p);
   static void deleteArray_VectorCompareResult(void *p);
   static void destruct_VectorCompareResult(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VectorCompareResult*)
   {
      ::VectorCompareResult *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::VectorCompareResult));
      static ::ROOT::TGenericClassInfo 
         instance("VectorCompareResult", "cpv_joint_Multi_noSystem.C", 80,
                  typeid(::VectorCompareResult), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &VectorCompareResult_Dictionary, isa_proxy, 4,
                  sizeof(::VectorCompareResult) );
      instance.SetNew(&new_VectorCompareResult);
      instance.SetNewArray(&newArray_VectorCompareResult);
      instance.SetDelete(&delete_VectorCompareResult);
      instance.SetDeleteArray(&deleteArray_VectorCompareResult);
      instance.SetDestructor(&destruct_VectorCompareResult);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VectorCompareResult*)
   {
      return GenerateInitInstanceLocal((::VectorCompareResult*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VectorCompareResult*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *VectorCompareResult_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::VectorCompareResult*)0x0)->GetClass();
      VectorCompareResult_TClassManip(theClass);
   return theClass;
   }

   static void VectorCompareResult_TClassManip(TClass* theClass){
      theClass->CreateAttributeMap();
      TDictAttributeMap* attrMap( theClass->GetAttributeMap() );
      attrMap->AddProperty("file_name","/exp/dune/app/users/cnguyen/lblpwgtools/CAFAna/tute/./cpv_joint_Multi_noSystem.C");
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_VectorCompareResult(void *p) {
      return  p ? new(p) ::VectorCompareResult : new ::VectorCompareResult;
   }
   static void *newArray_VectorCompareResult(Long_t nElements, void *p) {
      return p ? new(p) ::VectorCompareResult[nElements] : new ::VectorCompareResult[nElements];
   }
   // Wrapper around operator delete
   static void delete_VectorCompareResult(void *p) {
      delete ((::VectorCompareResult*)p);
   }
   static void deleteArray_VectorCompareResult(void *p) {
      delete [] ((::VectorCompareResult*)p);
   }
   static void destruct_VectorCompareResult(void *p) {
      typedef ::VectorCompareResult current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::VectorCompareResult

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 214,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace {
  void TriggerDictionaryInitialization_cpv_joint_Multi_noSystem_C_ACLiC_dict_Impl() {
    static const char* headers[] = {
"./cpv_joint_Multi_noSystem.C",
0
    };
    static const char* includePaths[] = {
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e15-prof/include",
"/exp/dune/app/users/cnguyen/lblpwgtools/CAFAna/build/Linux/include",
"/cvmfs/larsoft.opensciencegrid.org/products/boost/v1_66_0a/Linux64bit+3.10-2.17-e15-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_2a/Linux64bit+3.10-2.17-e15-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/xrootd/v4_8_0b/Linux64bit+3.10-2.17-e15-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/mysql_client/v5_5_58a/Linux64bit+3.10-2.17-e15/include",
"/cvmfs/larsoft.opensciencegrid.org/products/postgresql/v9_6_6a/Linux64bit+3.10-2.17-p2714b/include",
"/cvmfs/larsoft.opensciencegrid.org/products/pythia/v6_4_28k/Linux64bit+3.10-2.17-gcc640-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/gsl/v2_4/Linux64bit+3.10-2.17-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/fftw/v3_3_6_pl2/Linux64bit+3.10-2.17-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/clhep/v2_3_4_6/Linux64bit+3.10-2.17-e15-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/gmp/v6_2_1/Linux64bit+3.10-2.17/include",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e15-prof/etc",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e15-prof/etc/cling",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e15-prof/include",
"/cvmfs/larsoft.opensciencegrid.org/products/tbb/v2018_2a/Linux64bit+3.10-2.17-e15-prof/include",
"/scratch/workspace/canvas-products/vdevelop/e15/SLF7/prof/build/tbb/v2018_2a/Linux64bit+3.10-2.17-e15-prof/include",
"/usr/include/freetype2",
"/cvmfs/larsoft.opensciencegrid.org/products/root/v6_12_06a/Linux64bit+3.10-2.17-e15-prof/include",
"/exp/dune/app/users/cnguyen/lblpwgtools/CAFAna/tute/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "cpv_joint_Multi_noSystem_C_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
struct __attribute__((annotate(R"ATTRDUMP(file_name@@@/exp/dune/app/users/cnguyen/lblpwgtools/CAFAna/tute/./cpv_joint_Multi_noSystem.C)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(pattern@@@*)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$./cpv_joint_Multi_noSystem.C")))  VectorCompareResult;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "cpv_joint_Multi_noSystem_C_ACLiC_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "./cpv_joint_Multi_noSystem.C"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"", payloadCode, "@",
"GetAsimovName", payloadCode, "@",
"VectorCompareResult", payloadCode, "@",
"compareVectorsAsSets", payloadCode, "@",
"cpv_joint", payloadCode, "@",
"cpv_joint_Multi_noSystem", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("cpv_joint_Multi_noSystem_C_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_cpv_joint_Multi_noSystem_C_ACLiC_dict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_cpv_joint_Multi_noSystem_C_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_cpv_joint_Multi_noSystem_C_ACLiC_dict() {
  TriggerDictionaryInitialization_cpv_joint_Multi_noSystem_C_ACLiC_dict_Impl();
}
