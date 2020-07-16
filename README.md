# Algovault
This is a collection of algorithms, data structures and templates for competitive programming. In a nutshell, these are the goals of this project:
  1. Useful in online competitions like [Codeforces](https://codeforces.com/), [Google Code Jam](https://codingcompetitions.withgoogle.com/codejam), etc.
  2. Should be relatively simple to use as a black-box without compromising speed.
  3. Example usage and sufficient documentation.

Most codes are written in C/C++, and some are written in Python. While most implementations are stress-tested and cross-checked against various problems in several online judges, nonetheless they are not guranteed to be flawless and work in all cases. For bugs, refactoring and improvements, feel free to file an issue or a pull request as contributions are always welcome.

## Motivation
The initial goal was to only upload my personal templates online, nothing more :)

Shortly I stumbled upon this awesome project called [KACTL](https://github.com/kth-competitive-programming/kactl). KACTL is the ICPC team reference document of the KTH Royal Insitute of Technology made open source. This is great and probably the best team notebook I've seen, but the ICPC imposes a 25 page limit on the notebook and as a result there are many algorithms and variations helpful in online competitions that had to be left out.

Then there exists projects like [cosmos](https://github.com/OpenGenus/cosmos). While this is a nice initiative, but it's not catered towards competitive programming. And there are too many implementations of the same algorithms which is confusing and very specific.

There are some other amazing library like the ones created by [indy256](https://github.com/indy256/codelibrary), [bqi343](https://github.com/bqi343/USACO), [AhmadElsagheer](https://github.com/AhmadElsagheer/Competitive-programming-library), [kmyk](https://github.com/kmyk/competitive-programming-library), [OmarBazaraa](https://github.com/OmarBazaraa/Competitive-Programming) but they do share some common problems. Some focuses on basic/intermediate stuff and do not contain many advanced algorithms. Some are too modularized and cannot be copy pasted easily, for example FFT template has include reference from several other headers and are more like a third-party library. Additonally it's not intuitive how to use them as a black-box and often usage, common gotchas and specifics like time complexity and performance are left out. Most of these work are intended to be used as a personal collection, not focusing much on reusability.

It is for the above stated reasons I felt something like this could be helpful to the community. Note, for Python users, also check out this super cool project - [PyRival](https://github.com/cheran-senthil/PyRival).

## Getting Started
It's not exactly rocket science. Just copy paste any implementation and enjoy. You only need to keep in mind two things:
1. Most implementations should work with C++11 and some might require C++14 or higher.
2. For Python, use Python 3. Most codes should be compatible with PyPy too.

## Coding Style
1. No unncesessary macros or templates.
   * Because they suck and takes away readability.
   * Easier for others to copy paste code and use it with their own template.
   * Often unnecessarily lengthens the implementation.
2. Sufficient documentation is included.
   * Notes about usage, specifics and complexity in the header.
   * Documentation regarding commonly used/exposed methods.
   * Example usage.
3. Code should be easy to copy/paste. That is, no includes from custom or third-party headers even if this means duplication.
4. Implementation should be relatively simple, clean and easy to follow.
5. Performance is important. Speed should not be compromised significantly.
6. Indentation is 4 spaces everywhere and adequate.

## Future Work
1. This is still a work in progress so I'll port more code from my template over the time.
2. Add stress tests, automated unit tests and integration tests.
3. Refactor and simplify old implementations.
4. Add practice problems.
5. Add a comprehensive index.

## License
The project is licensed under the [GNU General Public License v3.0](https://github.com/sgtlaugh/algovault/blob/master/LICENSE). To summarize, this means that you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.
