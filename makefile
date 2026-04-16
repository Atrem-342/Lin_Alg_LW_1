CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -g
LDLIBS = -lm

APP_SRC = main.c gauss.c lu.c clock_get_time.c generate_matrix.c analysis.c matrix.c
TEST_SRC = tests.c gauss.c lu.c clock_get_time.c generate_matrix.c analysis.c matrix.c

all: matrix_app matrix_tests

matrix_app: $(APP_SRC)
	$(CC) $(CFLAGS) $(APP_SRC) $(LDLIBS) -o matrix_app

matrix_tests: $(TEST_SRC)
	$(CC) $(CFLAGS) $(TEST_SRC) $(LDLIBS) -o matrix_tests

test: matrix_tests
	./matrix_tests

clean:
	rm -f matrix_app matrix_tests
