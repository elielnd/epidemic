CC = clang++
OBJS_COMMON= matrix.o 
SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)


all: page_rank without_vaccine with_vaccine with_infection_vector

page_rank: $(OBJS_COMMON) page_rank.o
	$(CC)  -o $@ $^ $(CPPFLAGS)

without_vaccine: $(OBJS_COMMON) without_vaccine.o
	$(CC)  -o $@ $^ $(CPPFLAGS)

with_vaccine: $(OBJS_COMMON) with_vaccine.o
	$(CC)  -o $@ $^ $(CPPFLAGS)

with_infection_vector: $(OBJS_COMMON)  with_infection_vector.o
	$(CC)  -o $@ $^ $(CPPFLAGS)

matrix.o: src/matrix.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)

page_rank.o: src/page_rank.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)

without_vaccine.o: src/without_vaccine.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)

with_vaccine.o: src/with_vaccine.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)

with_infection_vector.o: src/with_infection_vector.cpp
	$(CC)  -o $@ -c  $< $(CPPFLAGS)





clean :	
	rm -f $(OBJ) *.o page_rank without_vaccine














# CC = clang++
# CPPFLAGS = 
# OUTPUT = spread_epidemic
# SRC_DIR = src
# SRC = $(wildcard $(SRC_DIR)/*.cpp)
# OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)


# $(SRC_DIR)/%.o : $(SRC_DIR)/%.cpp
# 	$(CC)  -o $@ -c  $< $(CPPFLAGS)  

# $(OUTPUT) : $(OBJ)
# 	$(CC)  -o $@ $^ $(CPPFLAGS)


# clean :
# 	rm -f $(OBJ) $(OUTPUT)