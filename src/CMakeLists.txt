SET(functions_SOURCES
    functions/Function.cpp
    functions/FunctionTime.cpp
    )

SET(integrators_SOURCES
    integrators/CVode.cpp
    integrators/DormandPrince45.cpp
    integrators/ExplicitRungeKuttaBase.cpp
    integrators/IntegratorStats.cpp
    integrators/OdeIntegratorBase.cpp
    integrators/RungeKutta23.cpp
    integrators/RungeKuttaBase.cpp
    integrators/RungeKuttaChebyshev.cpp
    integrators/RungeKuttaMerson45.cpp
    )

set(io_SOURCES
    io/DataWriter.cpp
    io/NetCDFReader.cpp
    )

SET(odesystem_SOURCES
    odesystem/AllenCahnEquation.cpp
    odesystem/CahnHilliardEquation.cpp
    odesystem/DegenerateCahnHilliardEquation.cpp
    odesystem/LoretiMarchEquation.cpp
    odesystem/ExplicitOde.cpp
    odesystem/MolOdeSystem.cpp
    )

SET(geometry_SOURCES
    geometry/RectangularDomain.cpp
    )

SET(utils_SOURCES
    utils/ConditionalOStream.cpp
    utils/Exceptions.cpp
    utils/JobIdentifier.cpp
    utils/LogStream.cpp
    utils/ParameterHandler.cpp
    utils/ProgressDisplay.cpp
    utils/OdeityApplication.cpp
    utils/Timer.cpp
    utils/Utilities.cpp
    utils/Vector.cpp
    )

SET(odeity_SOURCES
    ${integrators_SOURCES}
    ${odesystem_SOURCES}
    ${utils_SOURCES}
    ${functions_SOURCES}
    ${io_SOURCES}
    ${geometry_SOURCES}
    )

INCLUDE_DIRECTORIES(
    ${PROJECT_SOURCE_DIR}/vendor/sundials-2.4.0/include
    ${PROJECT_BINARY_DIR}/vendor/sundials-2.4.0/include
    )

# Add the build target for the static Odeity library
ADD_LIBRARY(odeity STATIC ${odeity_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(odeity PROPERTIES OUTPUT_NAME odeity CLEAN_DIRECT_OUTPUT 1)

# ADD_SUBDIRECTORY(bin)

#INSTALL(TARGETS odeity_integrators_static DESTINATION lib)

