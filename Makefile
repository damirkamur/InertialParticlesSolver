.PHONY: run tests clean

run:
	time julia --project -t 4 --optimize=3 main.jl