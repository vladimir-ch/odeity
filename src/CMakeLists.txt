# Copyright (c) 2010 Vladimir Chalupecky
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


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

