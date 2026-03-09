# ---------------------------------------------------------------------------
# Install C header files (always, regardless of SKBUILD)
# ---------------------------------------------------------------------------
install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/atormsd.h
    ${CMAKE_CURRENT_SOURCE_DIR}/conformsd.h
    DESTINATION include/molalignlib
)