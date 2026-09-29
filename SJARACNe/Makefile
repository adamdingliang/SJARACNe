#############################################################################
# Makefile for building: sjaracne
# Template: app
#############################################################################

UNAME_S := $(shell uname -s)

# CC
ifeq ($(UNAME_S), Darwin)
	CC := clang++ -arch x86_64
else
	CC := g++
endif

# Folders
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin
TOOLDIR := ../tools

# Targets
EXECUTABLE := sjaracne.exe
TARGET := $(TARGETDIR)/$(EXECUTABLE)
NULL_EXECUTABLE := apmi_null_generator.exe
NULL_TARGET := $(TARGETDIR)/$(NULL_EXECUTABLE)
NULL_OBJECT := $(BUILDDIR)/tools/apmi_null_generator.o

# Final Paths
INSTALLBINDIR := /usr/local/bin

# Code Listss
SRCEXT := cpp
HEADEXT := h
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SOURCES:.$(SRCEXT)=.o))

# Folder Lists
INCDIRS := $(shell find src/* -name '*.h' -exec dirname {} \; | sort | uniq)
INCLIST := $(patsubst src/%, -I src/%, $(INCDIRS))
BUILDLIST := $(patsubst src/%, $(BUILDDIR)/%, $(INCDIRS))

# Shared Compiler Flags
CFLAGS := -c -pipe -g -Wall -W -fPIC
INC := -I src $(INCLIST) -I /usr/local/include
LIB := -L /usr/local/lib

ifeq ($(UNAME_S), Linux)
	CFLAGS += -O3
else
	CFLAGS += -std=c++11 -stdlib=libc++ -O3
endif

.PHONY: all clean install distclean

all: $(TARGET) $(NULL_TARGET)

$(TARGET): $(OBJECTS)
	@mkdir -p $(TARGETDIR)
	@echo "Linking..."
	@echo "		Linking $(TARGET)"; $(CC) $^ -o $(TARGET) $(LIB)

$(NULL_TARGET): $(NULL_OBJECT) $(BUILDDIR)/apmi.o
	@mkdir -p $(TARGETDIR)
	@echo "Linking..."
	@echo "		Linking $(NULL_TARGET)"; $(CC) $^ -o $(NULL_TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDLIST) $(BUILDDIR)
	@echo "Compiling $<..."; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

$(NULL_OBJECT): $(TOOLDIR)/apmi_null_generator.cpp $(SRCDIR)/apmi.h
	@mkdir -p $(BUILDDIR)/tools
	@echo "Compiling $<..."; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning $(TARGET) and $(NULL_TARGET)..."; $(RM) -r $(BUILDDIR)/* $(TARGET) $(NULL_TARGET)

install: $(TARGET) $(NULL_TARGET)
	@echo "Installing $(EXECUTABLE) and $(NULL_EXECUTABLE)..."; cp $(TARGET) $(NULL_TARGET) $(INSTALLBINDIR)

distclean:
	@echo "Removing $(EXECUTABLE) and $(NULL_EXECUTABLE)"; rm $(INSTALLBINDIR)/$(EXECUTABLE) $(INSTALLBINDIR)/$(NULL_EXECUTABLE)
