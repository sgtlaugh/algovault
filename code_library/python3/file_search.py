import os
import sys
import glob


def run_script(path, extensions, patterns):
    """
    Small script to search and list all files containing any of the patterns as substring
    Useful to find old codes containing specific terms
    Can also be done with grep/find or shell commands but this should work cross-platform as well

    :param path: file path to search (not recursive)
    :param extensions: file extensions to search
    :param patterns: list of patterns to search for
    :return:
    """
    os.chdir(path)

    matched_files = []
    for ext in extensions:
        files = glob.glob('*{}'.format(ext))
        print('{} files found with extension {}'.format(len(files), ext))

        for file in files:
            sys.stdin = open(file)
            file_content = '\n'.join(line for line in sys.stdin)

            if any(pattern in file_content for pattern in patterns):
                matched_files.append(file)

    print('\n{} file(s) found containing any of {}\n'.format(len(matched_files), patterns))
    for i, matched_file in enumerate(matched_files):
        print('{:<4} {}'.format(str(i + 1) + '.', matched_file))


def main():
    path = '/home/sgtlaugh/Desktop/codes'
    extensions = ['.c', '.cpp', '.py']
    patterns = ['fft', 'ntt']

    run_script(path, extensions, patterns)


if __name__ == '__main__':
    main()
