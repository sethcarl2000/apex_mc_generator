// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME simc_hrs_dict
#define R__NO_DEPRECATION

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

// Header files passed as explicit arguments
#include "/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src/HRSTrack_t.hh"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *HRSTrack_t_Dictionary();
   static void HRSTrack_t_TClassManip(TClass*);
   static void *new_HRSTrack_t(void *p = nullptr);
   static void *newArray_HRSTrack_t(Long_t size, void *p);
   static void delete_HRSTrack_t(void *p);
   static void deleteArray_HRSTrack_t(void *p);
   static void destruct_HRSTrack_t(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HRSTrack_t*)
   {
      ::HRSTrack_t *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HRSTrack_t));
      static ::ROOT::TGenericClassInfo 
         instance("HRSTrack_t", "HRSTrack_t.hh", 15,
                  typeid(::HRSTrack_t), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HRSTrack_t_Dictionary, isa_proxy, 4,
                  sizeof(::HRSTrack_t) );
      instance.SetNew(&new_HRSTrack_t);
      instance.SetNewArray(&newArray_HRSTrack_t);
      instance.SetDelete(&delete_HRSTrack_t);
      instance.SetDeleteArray(&deleteArray_HRSTrack_t);
      instance.SetDestructor(&destruct_HRSTrack_t);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HRSTrack_t*)
   {
      return GenerateInitInstanceLocal((::HRSTrack_t*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HRSTrack_t*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HRSTrack_t_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HRSTrack_t*)nullptr)->GetClass();
      HRSTrack_t_TClassManip(theClass);
   return theClass;
   }

   static void HRSTrack_t_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HRSTrack_tcLcLHRSCoord_t_Dictionary();
   static void HRSTrack_tcLcLHRSCoord_t_TClassManip(TClass*);
   static void *new_HRSTrack_tcLcLHRSCoord_t(void *p = nullptr);
   static void *newArray_HRSTrack_tcLcLHRSCoord_t(Long_t size, void *p);
   static void delete_HRSTrack_tcLcLHRSCoord_t(void *p);
   static void deleteArray_HRSTrack_tcLcLHRSCoord_t(void *p);
   static void destruct_HRSTrack_tcLcLHRSCoord_t(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HRSTrack_t::HRSCoord_t*)
   {
      ::HRSTrack_t::HRSCoord_t *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HRSTrack_t::HRSCoord_t));
      static ::ROOT::TGenericClassInfo 
         instance("HRSTrack_t::HRSCoord_t", "HRSTrack_t.hh", 18,
                  typeid(::HRSTrack_t::HRSCoord_t), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HRSTrack_tcLcLHRSCoord_t_Dictionary, isa_proxy, 4,
                  sizeof(::HRSTrack_t::HRSCoord_t) );
      instance.SetNew(&new_HRSTrack_tcLcLHRSCoord_t);
      instance.SetNewArray(&newArray_HRSTrack_tcLcLHRSCoord_t);
      instance.SetDelete(&delete_HRSTrack_tcLcLHRSCoord_t);
      instance.SetDeleteArray(&deleteArray_HRSTrack_tcLcLHRSCoord_t);
      instance.SetDestructor(&destruct_HRSTrack_tcLcLHRSCoord_t);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HRSTrack_t::HRSCoord_t*)
   {
      return GenerateInitInstanceLocal((::HRSTrack_t::HRSCoord_t*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HRSTrack_t::HRSCoord_t*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HRSTrack_tcLcLHRSCoord_t_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HRSTrack_t::HRSCoord_t*)nullptr)->GetClass();
      HRSTrack_tcLcLHRSCoord_t_TClassManip(theClass);
   return theClass;
   }

   static void HRSTrack_tcLcLHRSCoord_t_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_HRSTrack_t(void *p) {
      return  p ? new(p) ::HRSTrack_t : new ::HRSTrack_t;
   }
   static void *newArray_HRSTrack_t(Long_t nElements, void *p) {
      return p ? new(p) ::HRSTrack_t[nElements] : new ::HRSTrack_t[nElements];
   }
   // Wrapper around operator delete
   static void delete_HRSTrack_t(void *p) {
      delete ((::HRSTrack_t*)p);
   }
   static void deleteArray_HRSTrack_t(void *p) {
      delete [] ((::HRSTrack_t*)p);
   }
   static void destruct_HRSTrack_t(void *p) {
      typedef ::HRSTrack_t current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HRSTrack_t

namespace ROOT {
   // Wrappers around operator new
   static void *new_HRSTrack_tcLcLHRSCoord_t(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::HRSTrack_t::HRSCoord_t : new ::HRSTrack_t::HRSCoord_t;
   }
   static void *newArray_HRSTrack_tcLcLHRSCoord_t(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::HRSTrack_t::HRSCoord_t[nElements] : new ::HRSTrack_t::HRSCoord_t[nElements];
   }
   // Wrapper around operator delete
   static void delete_HRSTrack_tcLcLHRSCoord_t(void *p) {
      delete ((::HRSTrack_t::HRSCoord_t*)p);
   }
   static void deleteArray_HRSTrack_tcLcLHRSCoord_t(void *p) {
      delete [] ((::HRSTrack_t::HRSCoord_t*)p);
   }
   static void destruct_HRSTrack_tcLcLHRSCoord_t(void *p) {
      typedef ::HRSTrack_t::HRSCoord_t current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HRSTrack_t::HRSCoord_t

namespace {
  void TriggerDictionaryInitialization_libsimc_hrs_Impl() {
    static const char* headers[] = {
"/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src/HRSTrack_t.hh",
nullptr
    };
    static const char* includePaths[] = {
"/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs",
"/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src",
"/u/group/halla/apps/ROOT/6.26-10/el9/RelWithDebInfo/include",
"/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src",
"/u/group/halla/apps/ROOT/6.26-10/el9/RelWithDebInfo/include",
"/u/group/halla/apps/ROOT/6.26-10/el9/RelWithDebInfo/include/",
"/w/halla-scshelf2102/apex/disk1/sethhall/apex_mc_generator/simc_hrs/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libsimc_hrs dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src/HRSTrack_t.hh")))  HRSTrack_t;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libsimc_hrs dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/work/halla/apex/disk1/sethhall/apex_mc_generator/simc_hrs/src/HRSTrack_t.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HRSTrack_t", payloadCode, "@",
"HRSTrack_t::HRSCoord_t", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libsimc_hrs",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libsimc_hrs_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libsimc_hrs_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libsimc_hrs() {
  TriggerDictionaryInitialization_libsimc_hrs_Impl();
}
