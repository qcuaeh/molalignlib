# ---------------------------------------------------------------------------
# Install C header files (always, regardless of SKBUILD)
# ---------------------------------------------------------------------------
install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/molalign.h
    ${CMAKE_CURRENT_SOURCE_DIR}/error_codes.h
    DESTINATION include/molalignlib
)