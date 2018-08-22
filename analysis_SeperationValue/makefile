# Declare
#
# C FLAGS
CC := g++
CC_FLAGS =-std=c++11 -O2 #-ansi -Wall -g
INC := -Iinc -I/home/kgameil/Code/boost/boost_1_64_0/
LIB := -L. -lMinuit -lPID -lInterp -L/home/kgameil/Code/boost/boost_1_64_0/stage/lib -lboost_iostreams
#ROOT FLAGS
R_LDFLAGS :=     `root-config --ldflags`
R_LIBS    :=     `root-config --glibs`
R_CFLAGS  :=     `root-config --cflags`
R_ALL     :=     $(R_LADFLAGS) $(R_LIBS) $(R_CFLAGS)

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
	make clean
	make
	"g++" -fPIC $(CC_FLAGS)$ $(R_ALL) $(LIB) $(INC)  $(OBJECTS) -shared -o libPID.so
