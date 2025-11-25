# 定义变量（补充必要编译选项，解决语法/兼容性问题）
CXX = g++
CXXFLAGS = -std=c++14
SRC_DIR = code
BIN_DIR = bin
TARGET = $(BIN_DIR)/VQBG  

# 列出所有源文件（保持不变）
SRCS = $(SRC_DIR)/common.cpp \
       $(SRC_DIR)/kmer_hash.cpp \
       $(SRC_DIR)/main.cpp \
       $(SRC_DIR)/sequence_graph.cpp \
       $(SRC_DIR)/utility.cpp \

# 列出所有对象文件（保持不变：code/xxx.o）
OBJS = $(SRCS:.cpp=.o)

# 头文件搜索路径（关键：若源文件包含 code/ 下的头文件，必须加这行）
INCLUDES = -I$(SRC_DIR) 


# 默认目标
all: $(TARGET)

# 规则：链接所有 .o 文件生成可执行文件（保持不变）
$(TARGET): $(OBJS)
	@mkdir -p $(BIN_DIR) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)  

# 规则：编译 .o 文件（核心修复！目标是 code/xxx.o，依赖 code/xxx.cpp）
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@  

# 清理目标（保持不变，清理 .o 和可执行文件）
clean:
	rm -f $(OBJS)
	rm -rf $(BIN_DIR) 

.PHONY: all clean