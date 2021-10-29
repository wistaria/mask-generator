file(GLOB PY_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.py")
foreach(name ${PY_SOURCES})
  configure_file(${name} ${name} COPYONLY)
endforeach(name ${PY_SOURCES})

pybind11_add_module(engine engine.cpp)
