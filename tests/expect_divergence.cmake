# A diverged run must stop with a non-zero exit and name the cycle. The test is
# on the process, because the old failure was a clean exit printing NaN.
execute_process(
    COMMAND ${PROGRAM} --voltage nan --cycles ${CYCLES} ${EXTRA}
    RESULT_VARIABLE code
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err)
if(code EQUAL 0)
    message(FATAL_ERROR "exited 0 on a NaN run:\n${out}")
endif()
if(NOT err MATCHES "solution diverged in cycle [0-9]+")
    message(FATAL_ERROR "exit ${code}, but no message naming the cycle:\n${err}")
endif()
