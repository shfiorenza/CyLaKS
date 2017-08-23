For program:
 .yaml files contain parameters; can be created using any text-editor
 .config files contain initial MT configs; created using ./buildmts
 .file files are the output of ./sim

For makefile:
# DEBUG mode
    make sim

# RELEASE mode (LIKE 100X FASTER COMPUTATIONAL TIME BRO!!!)
	make CFG=release sim
