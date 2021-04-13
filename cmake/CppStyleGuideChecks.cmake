# -----------------------------------------------------------------------------
#
# Copyright (C) 2021 CERN & Newcastle University for the benefit of the
# BioDynaMo collaboration. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------

# define scirpt to obtain source file lists
set(GET_SRC_FILES "${BUILD_SUPPORT_DIR}/get-src-files.sh")

# create helper target to check header files with clang-tidy
set(clang_tidy_header_helper ${CMAKE_BINARY_DIR}/clang_tidy_header_helper.cc)
file(WRITE ${clang_tidy_header_helper} "" ) # creates file
add_executable(clang-tidy-header-helper EXCLUDE_FROM_ALL ${clang_tidy_header_helper})
target_include_directories(clang-tidy-header-helper PUBLIC "${CMAKE_SOURCE_DIR}/test")

# -------------------- "make format*", "make show-format*" and "make check-format*" targets ------------
# @param file_filter: Filter that is passed to get-src-files scirpt (all|staged|changed)
function(add_clang_format_target make_target_id file_filter)
  if (${CLANG_FORMAT_FOUND})
    add_custom_target(${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-format.sh ${CMAKE_CURRENT_SOURCE_DIR} 
      ${CLANG_FORMAT_BIN} 1 `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-format on selected files and update them in-place")

    add_custom_target(show-${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-format.sh ${CMAKE_CURRENT_SOURCE_DIR} 
      ${CLANG_FORMAT_BIN} 2 `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-format on selected files and display differences")

    add_custom_target(check-${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-format.sh ${CMAKE_CURRENT_SOURCE_DIR} 
      ${CLANG_FORMAT_BIN} 0 `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-format on selected files. Fails if any file needs to be reformatted")
  endif()
endfunction(add_clang_format_target)

# -------------------- "make tidy*", "make show-tidy*" and "make check-tidy*" targets ------------
# @param file_filter: Filter that is passed to get-src-files scirpt (all|staged|changed)
function(add_clang_tidy_target make_target_id file_filter)
  if (${CLANG_TIDY_FOUND})
    add_custom_target(${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-tidy.sh ${CLANG_TIDY_BIN} ${clang_tidy_header_helper}  ${CMAKE_BINARY_DIR}/compile_commands.json 1
              `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h | grep -v -F -f ${PROJECT_SOURCE_DIR}/.clang-tidy-ignore`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-tidy on selected files and attempt to fix any warning automatically")

    add_custom_target(show-${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-tidy.sh ${CLANG_TIDY_BIN} ${clang_tidy_header_helper} ${CMAKE_BINARY_DIR}/compile_commands.json 2
              `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h | grep -v -F -f ${PROJECT_SOURCE_DIR}/.clang-tidy-ignore`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-tidy on selected files and display errors.")

    add_custom_target(check-${make_target_id}
      COMMAND ${BUILD_SUPPORT_DIR}/run-clang-tidy.sh ${CLANG_TIDY_BIN} ${clang_tidy_header_helper} ${CMAKE_BINARY_DIR}/compile_commands.json 0
              `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h | grep -v -F -f ${PROJECT_SOURCE_DIR}/.clang-tidy-ignore`
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      COMMENT "Run clang-tidy on selected files. Fails if errors are found.")
  endif()
endfunction(add_clang_tidy_target)

# -------------------- "make yapf*", "make show-yapf*" and "make check-yapf*" targets ------------
# @param file_filter: Filter that is passed to get-src-files scirpt (all|staged|changed)
function(add_yapf_target make_target_id file_filter)
  add_custom_target(${make_target_id}
    COMMAND ${BUILD_SUPPORT_DIR}/run-yapf.sh ${PROJECT_SOURCE_DIR} 1
            `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} py`
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run yapf on selected files and attempt to fix any warning automatically")

  add_custom_target(show-${make_target_id}
    COMMAND ${BUILD_SUPPORT_DIR}/run-yapf.sh ${PROJECT_SOURCE_DIR} 2
            `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} py`
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run yapf on selected files and display errors.")

  add_custom_target(check-${make_target_id}
    COMMAND ${BUILD_SUPPORT_DIR}/run-yapf.sh ${PROJECT_SOURCE_DIR} 0
            `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} py`
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run yapf on selected files. Fails if errors are found.")
endfunction(add_yapf_target)

# ------------------------------------------------------------------------------------------------

function(add_cpplint_target make_target_id file_filter)
  add_custom_target(${make_target_id}
    COMMAND ${BUILD_SUPPORT_DIR}/run-cpplint.sh `${GET_SRC_FILES} ${PROJECT_SOURCE_DIR} ${file_filter} cc h`
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Run cpplint on selected files. Fails if errors are found.")
endfunction(add_cpplint_target)

if (GIT_FOUND)
  add_custom_target(fetch-master
    COMMAND git fetch origin master
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Fetch latest changes from origin master. If you forked the project,
       make sure that it is synchronized with the biodynamo repository!")

  # Check if we have an origin/master branch
  execute_process(
  COMMAND bash -c "HAVE_ORIGIN=$(git branch -a | grep -c 'origin/master'); if [ $HAVE_ORIGIN -gt 0 ]; then echo 1; else echo 0; fi;"
    OUTPUT_VARIABLE HAVE_ORIGIN_REPOSITORY
  )

  # If origin/master is not set, then we disable some targets, more specifically:
  # - make format
  # - make tidy
  # - make check-cpplint-all
  if ("${HAVE_ORIGIN_REPOSITORY}" EQUAL "1")
    add_clang_format_target(format "changed")
    add_clang_tidy_target(tidy "changed")
    add_yapf_target(yapf "changed")
    add_cpplint_target(check-cpplint "changed")
  else()
    MESSAGE(WARNING "\nWe did not detect any remote origin/master branch. Therefore, targets \
like 'make format', 'make tidy' and 'make check-cpplint' will not be available. Use 'make format-all', 'make tidy-all' and 'make check-cpplint-all' instead. \
If you wish to use 'make format'/'make tidy'/'make check-cpplint' try to add an origin/master branch using 'git remote add'.\n")
  endif()

  add_clang_format_target(format-staged "staged")
  add_clang_tidy_target(tidy-staged "staged")
  add_yapf_target(yapf-staged "staged")
  add_cpplint_target(check-cpplint-staged "staged")

  # check submission
  add_custom_target(check-submission
    COMMAND "${BUILD_SUPPORT_DIR}/check-submission.sh" "${CMAKE_BINARY_DIR}"
    COMMENT "check-submission will build, run all tests, check formatting, code style, and generate documentation and coverage report")

  # fix submission
  add_custom_target(fix-submission
    COMMAND "${BUILD_SUPPORT_DIR}/fix-submission.sh" "${CMAKE_BINARY_DIR}"
    COMMENT "fix-submission will attempt to fix the reported issues using clang-format and clang-tidy.
             Failing tests and issues from cpplint must be fixed manually.")
endif()

add_clang_format_target(format-all "all")
add_clang_tidy_target(tidy-all "all")
add_yapf_target(yapf-all "all")
add_cpplint_target(check-cpplint-all "all")