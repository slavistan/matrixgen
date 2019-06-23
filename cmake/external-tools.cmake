# External tools. Run via `make <command>` from a build directory.

# clang-tidy
find_program(RUN_CLANG_TIDY_EXECUTABLE
             NAMES "run-clang-tidy.py"
             DOC "Path to the run-clang-tidy.py executable.")
if (RUN_CLANG_TIDY_EXECUTABLE)
  message(STATUS "run-clang-tidy found at '${RUN_CLANG_TIDY_EXECUTABLE}'.")
  add_custom_target(clang-tidy run-clang-tidy.py -checks="*" -header-filter=".*")
endif()

# 2. cppcheck
find_program(CPPCHECK_EXECUTABLE
             NAMES "cppcheck"
             DOC "Path to the cppcheck executable.")
if (CPPCHECK_EXECUTABLE)
  message(STATUS "cppcheck found at '${CPPCHECK_EXECUTABLE}'.")
  add_custom_target(cppcheck cppcheck --max-ctu-depth=3 --enable="all" --inconclusive --project=compile_commands.json)
endif()
