# Warning flags for our own targets only; dependencies are exempt (SPEC §6).

option(NTN_WARNINGS_AS_ERRORS "Treat warnings as errors (/WX) on our targets" OFF)

function(ntn_target_warnings target)
  target_compile_options(${target} PRIVATE /W4 /permissive- /utf-8 /EHsc)
  if(NTN_WARNINGS_AS_ERRORS)
    target_compile_options(${target} PRIVATE /WX)
  endif()
endfunction()
