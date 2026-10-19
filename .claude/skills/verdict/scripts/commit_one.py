import re
import subprocess
import sys

# commit_one.py <readme> <section header> <entry line> <solution path in repo> <message> [insert]
# Adds the ✅ link to the entry inside the section (or inserts the entry in natural order with "insert"),
# applies that one change to both HEAD's README (staged through the index only) and the working tree README,
# stages the solution file, and commits only those two paths, so another session's README edits stay out
readme, header, entry, solution, message = sys.argv[1:6]
insert = len(sys.argv) > 6 and sys.argv[6] == "insert"
link = f' [✅](../../{solution} "Accepted solution")'


def git(*args, inp=None):
    return subprocess.run(["git", "-C", repo, *args], input=inp, capture_output=True, text=True, check=True).stdout


def natural(line):
    m = re.match(r"- \[(.+?) - ", line)
    key = m.group(1) if m else line
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", key)]


def apply(text):
    lines = text.split("\n")
    start = lines.index(header)
    end = next((k for k in range(start + 1, len(lines)) if lines[k].startswith("#### ")), len(lines))
    body = [k for k in range(start + 1, end) if lines[k].startswith("- [")]
    if insert:
        assert not any(lines[k].startswith(entry) for k in body), "entry already present"
        pos = next((k for k in body if natural(lines[k]) > natural(entry)), body[-1] + 1 if body else start + 2)
        lines.insert(pos, entry + link)
    else:
        hits = [k for k in body if lines[k] == entry]
        assert len(hits) == 1, f"entry found {len(hits)} times without a link"
        lines[hits[0]] = entry + link
    return "\n".join(lines)


repo = subprocess.run(["git", "rev-parse", "--show-toplevel"], capture_output=True, text=True, check=True).stdout.strip()
new_head = apply(git("show", f"HEAD:{readme}"))
work_path = f"{repo}/{readme}"
with open(work_path) as f:
    updated = apply(f.read())
# Both versions are computed before anything is written, so a bad argument changes nothing
with open(work_path, "w") as f:
    f.write(updated)

blob = git("hash-object", "-w", "--stdin", inp=new_head).strip()
mode = git("ls-files", "-s", readme).split()[0]
git("update-index", "--cacheinfo", f"{mode},{blob},{readme}")
git("add", solution)
staged = set(git("diff", "--cached", "--name-only").split())
assert staged == {readme, solution}, f"unexpected staged paths: {sorted(staged)}"
print(git("diff", "--cached", "--stat").strip())
git("commit", "-q", "-m", message)
print(git("log", "--oneline", "-1").strip())
