# Algovault
This is a collection of algorithms, data structures and templates for competitive programming. In a nutshell, these are the goals of this project:
  1. Useful in online competitions like [Codeforces](https://codeforces.com/), [Google Code Jam](https://codingcompetitions.withgoogle.com/codejam), etc.
  2. Should be relatively simple to be used as a black-box without compromising speed.
  3. Example usage and sufficient documentation.

Most codes are written in C/C++, and some are written in Python. While most implementations are stress-tested and cross-checked against various problems in several online judges, nonetheless they are not guranteed to be flawless and work in all cases. For bugs, refactoring and improvements, feel free to file an issue or a pull request as contributions are always welcome.

## Motivation
The initial goal was to only upload my personal templates online, nothing more :-)

Then I stumbled upon this awesome project called [KACTL](https://github.com/kth-competitive-programming/kactl). KACTL is the ICPC team reference document of the KTH Royal Insitute of Technology made open source. This is great and probably the best team notebook I've seen, but the ICPC have a 25 page limit on the notebook and as a result there are many algorithms extremely useful in online competitions that had to be left out.

Then there exists projects like [cosmos](https://github.com/OpenGenus/cosmos). While this is a nice initiative, but it's not catered towards competitive programming. And there are too many implementations of the same algorithms which is confusing and very specific.

There are some amazing projects like the ones created by [bqi343](https://github.com/bqi343/USACO), [AhmadElsagheer](https://github.com/AhmadElsagheer/Competitive-programming-library), [kmyk](https://github.com/kmyk/competitive-programming-library), [OmarBazaraa](https://github.com/OmarBazaraa/Competitive-Programming) but they do share some common problems. Some focuses on basic/intermediate stuff and do not contain many advanced algorithms. Some are too modularized and cannot be copy pasted easily, for example FFT template has include reference from several other headers and are more like a Library than templates. Moreover it's not intuitive how to use them as a black-box and often usage, common gotchas and specifics like time complexity and performance are left out.

 
