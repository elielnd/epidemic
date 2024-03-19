CC = clang++
CPPFLAGS = 
LDFLAGS =-fopenmp
OUTPUT = page_rank
SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)


$(SRC_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) -o $@ -c $< $(CPPFLAGS) $(LDFLAGS)

$(OUTPUT) : $(OBJ)
	$(CC) -o $@ $^ $(CPPFLAGS) $(LDFLAGS)

clean :
	rm -f $(OBJ) $(OUTPUT)