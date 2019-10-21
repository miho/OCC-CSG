IF(UNIX)
  set(OS_SPECIFIC "pthread")
ENDIF(UNIX)  

CreateExecutable(occ-csg.cpp

#${OpenCASCADE_LIBRARIES}

TKSTL       # stl file support

# order of the library might matter for static linking
# remember to put base libs such as kernel and math at the end

TKXDESTEP
TKBinXCAF
TKXmlXCAF
TKXDEIGES
TKXCAF
TKIGES
TKSTEP
TKSTEP209
TKSTEPAttr
TKSTEPBase
TKXSBase 
TKStd 
TKStdL
TKXml 
TKBin 
TKXmlL 
TKBinL 
TKCAF 
TKXCAF 
TKLCAF 
TKCDF 
TKMeshVS 
TKOpenGl 
TKV3d
TKService
TKXMesh 
TKMesh 
TKOffset 
TKFeat
TKFillet 
TKHLR 
TKBool
TKBO
TKShHealing
TKPrim
TKTopAlgo
TKGeomAlgo
TKBRep
TKGeomBase
TKG3d
TKG2d
TKMath
TKernel

${FreeType}

${OS_SPECIFIC}

)

IF(UNIX AND NOT APPLE)
  link_target_libraries(occ-csg ${FontConfig})
ENDIF()  

install (TARGETS occ-csg DESTINATION bin)



#TKPrim     # primitives, such as sphere, cone etc.
#TKBool     # boolean operators
#TKSTEP     # step file support
#TKV3d      # font support
#TKFillet   # rounded edges (fillets)
#TKGeomAlgo # 
#TKMath     # math support
#TKernel    # base
