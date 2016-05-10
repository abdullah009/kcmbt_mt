
CC			= g++
SRC_DIR		= src
EXEC_DIR	= bin
SRC			= main_kcmbt.cc
EXEC		= kcmbt
DUMP_SRC	= decode.cc
DUMP_EXEC	= kcmbt_dump
QUERY_SRC	= sort_kmer.cc
QUERY_EXEC	= kcmbt_query
CFLAGS		= -Wno-unused-result -std=c++11 -O3 -pthread
CLINK		= -lz 

all: kcmbt kcmbt_dump kcmbt_query

kcmbt:
	mkdir -p bin
	$(CC) $(CFLAGS) $(SRC_DIR)/$(SRC) -o $(EXEC_DIR)/$(EXEC) $(CLINK)
	
kcmbt_dump:
	$(CC) $(CFLAGS) $(SRC_DIR)/$(DUMP_SRC) -o $(EXEC_DIR)/$(DUMP_EXEC) $(CLINK)
	
kcmbt_query:
	$(CC) $(CFLAGS) $(SRC_DIR)/$(QUERY_SRC) -o $(EXEC_DIR)/$(QUERY_EXEC)

clean:
	rm -f $(EXEC_DIR)/$(EXEC) $(EXEC_DIR)/$(DUMP_EXEC) $(EXEC_DIR)/$(QUERY_EXEC) *out* kmer_list.txt


