# Git & GitHub from the Command Line — Workshop

## Overview
By the end of this 2-hour workshop, every student has a GitHub repository for their MKR WiFi sensor project, has pushed at least 5 meaningful commits to it, and can do each step from the terminal.

|  |  |
| --- | --- 
| Prerequisites | VS Code installed; an Arduino sketch that compiles (a blink or WiFi-scan sketch is fine) |
| Bring | Laptop, MKR, Nicla|
| Output | A repo URL submitted on Canvas, containing code, README and a history of commits |

**How the workshop runs.** Every task is done with git commands typed into a terminal. The VS Code integrated terminal (opened with Ctrl+\` / Cmd+\`) is the most convenient, but any terminal works.

## Part 0 — Setup (15 min)

Everyone needs git installed, a GitHub account, and a way for git to authenticate with GitHub before anything else works. You've hopefully done this already as part of Workshop 0, but if you haven't already, the instructions are below:gi

1. **Install git.** Windows: [git-scm.com](https://git-scm.com/downloads), accept the defaults (this also installs Git Credential Manager). macOS: run `git --version` in a terminal and accept the Xcode tools prompt. Check with `git --version`.
1. **Tell git who you are** (CLI, once per computer). Use the same email as your GitHub account:
```bash
git config --global user.name "Ada Lovelace"
git config --global user.email "ada@liverpool.ac.uk"
git config --global init.defaultBranch main
```
1. **Authenticate with GitHub.** GitHub no longer accepts account passwords over the command line. When `git push` first asks for credentials, either let the credential-manager browser window appear and approve it, or install the [GitHub CLI](https://cli.github.com) and run `gh auth login`.

**Checkpoint:** `git config --global --list` shows your name and email, and (if you installed it) `gh auth status` shows you are logged in.

## Part 1 — Five ideas you need (10 min)
Git is a lab notebook for code: every commit is a dated, signed entry you can always return to; GitHub is the shared copy of that notebook online.

1. **Repository (repo):** a project folder that git is tracking. The history lives in a hidden `.git` folder inside it.
2. **Working → staged → committed:** you edit files (working), choose which changes go into the next entry (staged), then record them with a message (commit). `git status` shows which files are in each state.
3. **Commit:** Changes to the code that you're adding, plus a message saying *why*. Make small, frequent commits rather than one giant commit the night before the deadline. You're assessed on your commit record - both as evidence of incremental improvement and of teamwork.
4. **Remote:** the copy on GitHub, called `origin`. **Push** sends your commits up. **Pull** brings others' (or your other laptop's) commits down.
5. **Branch:** a parallel line of work. `main` is the version that works; try risky changes on a separate branch.

The flow you should be repeating during your whole project: **edit → stage → commit → push**.

## Part 2 — Put your project on GitHub (20 min)

You will turn your existing sketch folder into a repo and publish it as a **private** GitHub repository.

**Task 2.1 — Set up the folder.** Open a terminal in your project folder: in VS Code use *File → Open Folder* and the integrated terminal will start there, or `cd` into the folder from any terminal. Check with `pwd`.

**Task 2.2 — Initialise and publish.**

| Step | Command |
| --- | --- |
| Start tracking | `git init` |
| Stage everything | `git add .` |
| First commit | `git commit -m "Initial commit: WiFi sensor sketch"` |
| Create GitHub repo and push | Create an empty private repo on github.com (no README), then `git remote add origin https://github.com/<user>/<repo>.git` and `git push -u origin main` |

**Checkpoint:** refresh your repo page on github.com and see your sketch and `.gitignore`.

## Part 3 — The everyday loop (25 min)

Most of your git use is one loop: edit, look at what changed, stage, commit, push.

As well as tracking what we do, we can also give git a list of files that we want it to ignore - this is important if you have files containing personal information (e.g., WiFi passwords).

**A `.gitignore` for this project.** Create a file called `.gitignore` in the project root. Files listed here are never tracked:

```
# WiFi credentials — never commit
arduino_secrets.h

# OS and editor clutter
.DS_Store
Thumbs.db
.vscode/c_cpp_properties.json

# Python analysis
__pycache__/
.ipynb_checkpoints/
.venv/
```

**Task 3.1 — The loop.** Change the sensor read interval in your sketch, then:

```bash
git status                 # what has changed?
git diff                   # show the changes line by line
git add sensor_sketch.ino  # stage one file (or: git add .)
git commit -m "Change sensor read interval"
git push
git log --oneline          # see your history
```

**Task 3.2 — Again, from memory.** Add a comment block at the top of the sketch describing the board and sensor, and run through the same loop without looking at the commands above.

**Good commit messages** say what and why in under \~60 characters, written as a command:

| Weak | Better |
| --- | --- |
| `stuff` | `Add MQTT publish every 30 s` |
| `fixed it` | `Fix WiFi reconnect after router drop` |
| `final final v2` | `Calibrate temperature offset (−1.2 °C)` |

Aim for a commit each time something works or each lab session, whichever comes first.

## Part 4 — Branches, pulling, and undoing (25 min)

Branches let you try a change without breaking the version that works; pulling keeps two machines in step; a handful of commands undo most mistakes.

**Task 4.1 — Try a feature on a branch.** Add a second sensor reading (or change the upload format) on a new branch, merge it back once it works.

| Step | Command |
| --- | --- |
| New branch | `git switch -c add-humidity` |
| Commit and push work | `git add .` → `git commit -m "..."` → `git push -u origin add-humidity` |
| Merge into main | `git switch main` → `git merge add-humidity` → `git push` |

Alternatively, open a **pull request** on github.com and merge it there; this leaves a visible record of the change for your marker.

**Task 4.2 — Pull a change made elsewhere.** On github.com, edit `README.md` with the pencil icon and commit. Back in your terminal, run `git pull`. Always pull before you start work on a second computer.

**Task 4.3 — Undo.**

| Situation | Command |
| --- | --- |
| Throw away uncommitted edits to a file | `git restore <file>` |
| Unstage a file | `git restore --staged <file>` |
| Fix the last commit message (not yet pushed) | `git commit --amend -m "Better message"` |
| Undo a pushed commit safely | Find the commit with `git log --oneline`, then `git revert <commit-id>` |

**Merge conflicts.** If you and GitHub both changed the same lines, `git pull` stops and reports a conflict. `git status` lists the affected files; inside each, git marks the clash between `<<<<<<<`, `=======` and `>>>>>>>` lines. Edit the file to keep the version you want, delete the markers, then `git add` and `git commit`. Pulling before you start each session avoids most conflicts.

**README checklist** (write it in Markdown; VS Code previews with Ctrl+Shift+V):

- [ ] One-paragraph summary: what it measures, and why
- [ ] Hardware list and a wiring diagram or photo
- [ ] Libraries needed, with versions
- [ ] How to build and upload, including creating `arduino_secrets.h`
- [ ] Where the data goes (serial, cloud dashboard, MQTT) and a sample of it
- [ ] Known limitations and what you would do next

**Task 5.1:** Restructure your repo to the layout above, write a first README, commit and push. Check on github.com that it renders.

## Before you submit

- [ ] Repo is private and your module lead's GitHub account is added under *Settings → Collaborators*
- [ ] `git status` says *nothing to commit, working tree clean* and everything is pushed
- [ ] Code compiles from a fresh clone after copying the secrets template
- [ ] No passwords, keys or large build files anywhere in the repo
- [ ] README complete; history shows regular, descriptive commits across the project

## Cheat sheet

| Task | Command |
| --- | --- |
| Status | `git status` |
| See changes | `git diff` |
| Stage | `git add <file>` / `git add .` |
| Commit | `git commit -m "msg"` |
| Push / pull | `git push` / `git pull` |
| History | `git log --oneline` |
| New branch | `git switch -c <name>` |
| Clone | `git clone <url>` |

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| *Please tell me who you are* | Run the two `git config --global` lines in Part 0 |
| Push asks for a password, then fails | Passwords are not accepted: approve the credential-manager browser window, or run `gh auth login` |
| *rejected — fetch first* / *non-fast-forward* | Someone (or GitHub's web editor) changed the remote: `git pull`, resolve any conflict, then `git push` |
| *fatal: not a git repository* | Your terminal is in the wrong folder: `cd` into your project folder (or run `git init` if you haven't yet) |
| *src refspec main does not match any* | No commits yet: commit before your first push |
