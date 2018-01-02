Simulated Annealing Single-Cell (SASC) -- cancer progression inference
===================

Compile
--------

*OSX*
```bash
    gcc -o sasc *.c -std=c99
```

*Linux*
```bash
    gcc -o sasc *.c -std=gnu99 -g -lm
```

Run with:

```bash
    ./sasc -i data/gawad2.txt -m 16 -n 115 -a 0.3 -b 0.1 -k 2 -e data/gawad2_mut.txt
    ./sasc -i data/hou.txt -m 18 -n 58 -a 0.43 -b 0.0000604 -k 2 -e data/hou_mut.txt
```