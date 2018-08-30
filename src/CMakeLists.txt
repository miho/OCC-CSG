CreateExecutable(occ-csg.cpp

${OpenCASCADE_LIBRARIES}

#TKernel    # base
#TKSTL      # stl file support
#TKMath     # math support
#TKPrim     # primitives, such as sphere, cone etc.
#TKBool     # boolean operators
#TKSTEP     # step file support
#TKV3d      # font support
#TKFillet   # rounded edges (fillets)

        freetype

)

install (TARGETS occ-csg DESTINATION bin)
