import re
import sys

# api_order.py <template>...: prints each template whose struct defines the header's calls out of header order.
# Only member calls "obj.name(" count, grouped per obj: free functions may define helpers first, and a call defined
# outside the struct that holds most of its group (an inherited one) is left out
API_LINE = re.compile(r"^ \* ([a-z_]\w*)\.([A-Za-z_]\w*)(?:<[^>]*>)?\(")


def header_calls(lines):
    calls = []
    for line in lines:
        if line.startswith("***/"):
            break
        m = API_LINE.match(line)
        if m and m.groups() not in calls:
            calls.append(m.groups())
    return calls


def struct_of(lines):
    owner, current = [], -1
    for i, line in enumerate(lines):
        if re.match(r"^\s*struct \w+.*\{\s*$", line):
            current = i
        owner.append(current)
    return owner


def first_definition(lines, name):
    pattern = re.compile(r"^\s*(?:template\s*<[^>]*>\s*)?(?:static |inline |constexpr |friend )*[\w:<>,*& ]*?[\s*&]" + re.escape(name) + r"\s*\(")
    code = next((i + 1 for i, line in enumerate(lines) if line.startswith("***/")), 0)
    for i in range(code, len(lines)):
        line = lines[i]
        if "{" in line and pattern.match(line) and not line.strip().startswith(("return", "if", "for", "while")):
            return i
    return None


bad = 0
for path in sys.argv[1:]:
    lines = open(path).read().split("\n")
    owner = struct_of(lines)
    groups = {}
    for obj, name in header_calls(lines):
        groups.setdefault(obj, []).append(name)

    for obj, names in groups.items():
        defined = [(first_definition(lines, n), n) for n in names]
        defined = [(i, n) for i, n in defined if i is not None]
        if not defined:
            continue
        structs = [owner[i] for i, _ in defined]
        home = max(set(structs), key=structs.count)
        defined = [(i, n) for i, n in defined if owner[i] == home]
        if [n for _, n in defined] != [n for _, n in sorted(defined)]:
            bad += 1
            print(f"{path} ({obj}): header order {[n for _, n in defined]}, code order {[n for _, n in sorted(defined)]}")
sys.exit(1 if bad else 0)
