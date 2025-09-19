# 这里是目录喵～
DIR_OBJ = obj
DIR_SRC = src
DIR_INC = inc

# 要在 src 文件夹里找所有的源文件喵～
CPP_SRCS = $(wildcard $(DIR_SRC)/*.cpp)

# 把目标文件变成 $(patsubst 原模式，目标模式，文件列表) 的样子喵～
AOBJS=$(patsubst $(DIR_SRC)/%.cpp, $(DIR_OBJ)/%.o, $(CPP_SRCS))

# 头文件和依赖文件在这里喵～
DEPS=$(wildcard $(DIR_INC)/*.hpp)

# 要加上库和头文件的路径喵～
LIBS=
INCLUDES= -I./$(DIR_INC)

PROG = gmres

$(DIR_OBJ)/%.o:$(DIR_SRC)/%.cpp
	@mkdir -p $(DIR_OBJ)
	$(CXX) -c $(CXX_CFLAGS) $(INCLUDES)  $< -o $@

$(PROG): $(AOBJS)
	$(CXX) $(CXX_LDFLAGS) $(LIBS) $(INCLUDES)  $^ -o $@

.PHONY:clean
clean:
	rm -f $(DIR_OBJ)/*.o $(PROG)

.PHONY:format
format:
	clang-format -i $(DEPS) $(CPP_SRCS)

