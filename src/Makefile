# ADJUST THESE VARIABLES AS NEEDED TO SUIT YOUR COMPUTING ENVIRONMENT
# -------------------------------------------------------------------
# This specifies the C++ compiler used.
CXX = g++

# If the Boost and Eigen C++ include files are in non-standard
# locations, their locations should be declared here. For example, if
# you installed Boost and Eigen in the "MAPS" conda environment
# (following the instructions in the README), CXXFLAGS would most
# likely look like this:
#
#   CXXFLAGS = -I$(HOME)/anaconda/envs/MAPS/include -O3 -DNDEBUG
# 
CXXFLAGS = -O3 -DNDEBUG

# If the Boost C++ library is installed in a non-standard location,
# the location should be declared here. For example, if you installed
# Boost in the "MAPS" conda environment (following the instructions in
# the README), LDFLAGS would probably look like this:
#  
#   LDFLAGS = -L$(HOME)/anaconda/envs/MAPS/lib \
#             -lboost_system -lboost_program_options \
#             -lboost_filesystem
#
LDFLAGS = -lboost_system -lboost_program_options -lboost_filesystem

# AVOID EDITING ANYTHING BELOW THIS LINE
# --------------------------------------
# These two variables define the targets: the executable, and the
# object files.
OUTPUT = runeems2
OBJS = runeems2.o eems2.o util.o mcmc.o draw.o habitat.o graph.o

# These lines below define the rules for generating the targets.
all: $(OUTPUT)

$(OUTPUT): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OUTPUT) $(OBJS)

