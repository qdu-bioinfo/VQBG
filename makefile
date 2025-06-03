# 定义变量
CXX = g++
SRC_DIR = code
BIN_DIR = bin
TARGET = $(BIN_DIR)/VQBG

# 列出所有源文件
SRCS = $(SRC_DIR)/common.cpp \
       $(SRC_DIR)/kmer_hash.cpp \
       $(SRC_DIR)/main.cpp \
       $(SRC_DIR)/sequence_graph.cpp \
       $(SRC_DIR)/utility.cpp \

# 列出所有对象文件
OBJS = $(SRCS:.cpp=.o)

# 默认目标
all: $(TARGET)

# 规则来构建目标
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR) 
	$(CXX) -o $@ $(OBJS)

# 规则来构建每个对象文件
%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理目标
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
