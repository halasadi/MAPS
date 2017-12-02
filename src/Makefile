EIGEN_INC = /usr/local/include/eigen3
BOOST_LIB = /usr/local/opt/boost@1.57/lib
BOOST_INC = /usr/local/opt/boost@1.57/include

EXE = runeems2
OBJ = runeems2.o eems2.o util.o mcmc.o draw.o habitat.o graph.o
CXXFLAGS = -I${BOOST_INC} -I${EIGEN_INC} -O3 -DNDEBUG -O3
LDFLAGS = -L${BOOST_LIB} \
	-lboost_system \
	-lboost_program_options \
	-lboost_filesystem

$(EXE): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE) *~
