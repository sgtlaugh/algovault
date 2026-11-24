"""member_order.py <files>: prints file:line for each struct field or type declared after the struct's first method

Nested and indented structs are checked too. A method is a member line with a parameter list that opens a body,
on its own line or as a one-line body like int f(){ return 1; }.
"""
import re
import sys

STRUCT = re.compile(r'^( *)struct \w+( : [^{]+)?\{$')
METHOD = re.compile(r'^[^(=]*\(.*\{(.*\})?$')


def check(path):
    stack = []  # [indent of the struct's members, seen a method]
    for number, line in enumerate(open(path, encoding='utf-8'), 1):
        line = line.rstrip('\n')
        indent = len(line) - len(line.lstrip(' '))
        while stack and line.strip() and indent < stack[-1][0]:
            stack.pop()
        if stack and indent == stack[-1][0] and line[indent] not in '/}':
            body = line[indent:]
            if METHOD.match(body):
                stack[-1][1] = True
            elif stack[-1][1] and body.endswith(';') and '(' not in body and not re.match(r'(return|friend) ', body):
                print(f'{path}:{number}')
        match = STRUCT.match(line)
        if match:
            stack.append([len(match.group(1)) + 4, False])


for path in sys.argv[1:]:
    check(path)
