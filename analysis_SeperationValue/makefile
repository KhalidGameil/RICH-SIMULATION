# Declare
#

# C FLAGS
CC := g++
CC_FLAGS = -O2
INC = -Iinc -IGTEInclude
LIB = -L. -lPID -lMinuit -lInterp -lGTEngine -lprofiler

#ROOT FLAGS
R_LDFLAGS =     `root-config --ldflags`
R_LIBS    =     `root-config --glibs`
R_CFLAGS  =     `root-config --cflags`
R_ALL     =     $(R_LADFLAGS) $(R_LIBS) $(R_CFLAGS)

#BOOST FLAGS
B_FLAGS = #`-L/usr/local/include/boost`

# File names
EXEC = run
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

$(EXEC): $(OJECTS)
	$(CC)  $(CC_FLAGS) -fPIC $(LDFLAGS) $(INC) $(R_ALL) $(B_FLAGS) -I $(LIB) -c $(SOURCES)
	$(CC)  $(INC) $(P_CC_FLAG) $(CC_FLAGS) $(R_ALL) $(G_FLAGS) $(B_FLAGS) $(OBJECTS) $(LIB) -o final.out


# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

lib:
	"g++" -fPIC $(P_CC_FLAGS)$ $(R_ALL)  $(OBJECTS) -o libPID.dylib libInterp.dylib
