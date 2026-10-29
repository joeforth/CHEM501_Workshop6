# Git & GitHub in VS Code — Workshop

## Overview
Note this workshop borrows heavily from material [here](https://code.visualstudio.com/docs/sourcecontrol/quickstart).

By the end of this 2-hour workshop, every student has a GitHub repository for their MKR WiFi sensor project, has pushed at least 5 meaningful commits to it, and can do each step both from the VS Code Source Control panel and from the terminal.

|  |  |
| --- | --- 
| Prerequisites | VS Code installed; an Arduino sketch that compiles (a blink or WiFi-scan sketch is fine) |
| Bring | Laptop, MKR, Nicla|
| Output | A repo URL submitted on Canvas, containing code, README and a history of commits |

**How the workshop runs.** Every task is shown twice: a **GUI** route (VS Code clicks) and a **CLI** route (commands in the VS Code integrated terminal, opened with Ctrl+\` / Cmd+\`). Do each task with the GUI first, then repeat the next one with the CLI. Both drive the same git underneath, so you can mix them freely.

## Part 0 — Setup (15 min)

Everyone needs git installed, a GitHub account, and VS Code signed in to it before anything else works. You've hopefully done this already as part of Workshop 0, but if you haven't already, the instructions are below:gi

1. **Install git.** Windows: [git-scm.com](https://git-scm.com/downloads), accept the defaults (this also installs Git Credential Manager). macOS: run `git --version` in a terminal and accept the Xcode tools prompt. Check with `git --version`.
1. **Tell git who you are** (CLI, once per computer). Use the same email as your GitHub account:
%
%
```bash
git config --global user.name "Ada Lovelace"
git config --global user.email "ada@liverpool.ac.uk"
git config --global init.defaultBranch main
```
%
1. **Sign VS Code in to GitHub.** GUI: click the Accounts icon (person, bottom-left of the Activity Bar) → *Sign in with GitHub* → approve in the browser. This also lets VS Code push without asking for a password.
1. **Optional CLI login.** GitHub no longer accepts account passwords over the command line. If a terminal `git push` asks for credentials, either let the credential-manager browser window appear, or install the [GitHub CLI](https://cli.github.com) and run `gh auth login`.

**Checkpoint:** `git config --global --list` shows your name and email, and the Accounts icon shows your GitHub username.

## Part 1 — Five ideas you need (10 min)
Git is a lab notebook for code: every commit is a dated, signed entry you can always return to; GitHub is the shared copy of that notebook online.

1. **Repository (repo):** a project folder that git is tracking. The history lives in a hidden `.git` folder inside it.
2. **Working → staged → committed:** you edit files (working), choose which changes go into the next entry (staged), then record them with a message (commit). In VS Code, changed files appear under *Changes*; clicking **+** stages them.
3. **Commit:** Changes to the code that you're adding, plus a message saying *why*. Make small, frequent commits rather than one giant commit the night before the deadline. You're assessed on your commit record - both as evidence of incremental improvement and of teamwork.
4. **Remote:** the copy on GitHub, called `origin`. **Push** sends your commits up. **Pull** brings others' (or your other laptop's) commits down.
5. **Branch:** a parallel line of work. `main` is the version that works; try risky changes on a separate branch.

The flow you will repeat all term: **edit → stage → commit → push**.

## Part 2 — Put your project on GitHub (20 min)

You will turn your existing sketch folder into a repo and publish it as a **private** GitHub repository.

**Task 2.1 — Set up the folder.** Open your project folder in VS Code (*File → Open Folder*).

**Task 2.2 — Initialise and publish.**

| Step | GUI (VS Code) | CLI (integrated terminal) |
| --- | --- | --- |
| Start tracking | Source Control icon (Ctrl+Shift+G) → *Initialize Repository* | `git init` |
| Stage everything | Hover *Changes* → **+** | `git add .` |
| First commit | Type "Initial commit: WiFi sensor sketch" in the message box → **Commit** | `git commit -m "Initial commit: WiFi sensor sketch"` |
| Create GitHub repo and push | **Publish Branch** → choose *Publish to GitHub private repository* | Create an empty private repo on github.com (no README), then `git remote add origin https://github.com/<user>/<repo>.git` and `git push -u origin main` |

**Checkpoint:** refresh your repo page on github.com and see your sketch and `.gitignore`.

## Part 3 — The everyday loop (25 min)

Most of your git use is one loop: edit, look at what changed, stage, commit, push.

**A `.gitignore` for this project.** Create a file called `.gitignore` in the project root. Files listed here are never tracked:

```
# WiFi credentials — never commit
arduino_secrets.h

# Build output
build/
.pio/
*.bin
*.elf
*.hex

# OS and editor clutter
.DS_Store
Thumbs.db
.vscode/c_cpp_properties.json

# Python analysis
__pycache__/
.ipynb_checkpoints/
.venv/
```

**Task 3.1 — Loop with the GUI.** Change the sensor read interval in your sketch. In Source Control, click the file to see a side-by-side diff (red removed, green added). Stage it with **+**, write a message, **Commit**, then **Sync Changes** (or the ↑ arrow in the status bar).

**Task 3.2 — Loop with the CLI.** Add a comment block at the top of the sketch describing the board and sensor, then:

```bash
git status                 # what has changed?
git diff                   # show the changes line by line
git add sensor_sketch.ino  # stage one file (or: git add .)
git commit -m "Document board and sensor in header"
git push
git log --oneline          # see your history
```

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

| Step | GUI (VS Code) | CLI |
| --- | --- | --- |
| New branch | Click the branch name (`main`) in the status bar → *Create new branch* → `add-humidity` | `git switch -c add-humidity` |
| Commit and push work | As in Part 3; first push offers **Publish Branch** | `git add .` → `git commit -m "..."` → `git push -u origin add-humidity` |
| Merge into main | Switch to `main` via status bar → Command Palette → *Git: Merge Branch* → `add-humidity` | `git switch main` → `git merge add-humidity` → `git push` |

Alternatively, open a **pull request** on github.com and merge it there; this leaves a visible record of the change for your marker.

**Task 4.2 — Pull a change made elsewhere.** On github.com, edit `README.md` with the pencil icon and commit. Back in VS Code: **Sync Changes** (GUI) or `git pull` (CLI). Always pull before you start work on a second computer.

**Task 4.3 — Undo.**

| Situation | GUI | CLI |
| --- | --- | --- |
| Throw away uncommitted edits to a file | Right-click file under *Changes* → *Discard Changes* | `git restore <file>` |
| Unstage a file | **−** next to it under *Staged Changes* | `git restore --staged <file>` |
| Fix the last commit message (not yet pushed) | *Commit* dropdown (…) → *Commit Staged (Amend)* | `git commit --amend -m "Better message"` |
| Undo a pushed commit safely | Timeline view (Explorer) to inspect; revert via CLI | `git revert <commit-id>` |

**Merge conflicts.** If you and GitHub both changed the same lines, VS Code highlights the clash with *Accept Current / Accept Incoming / Accept Both* buttons. Pick one, save, stage and commit. Pulling before you start each session avoids most conflicts.

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

| Task | VS Code GUI | CLI |
| --- | --- | --- |
| Status | Source Control panel | `git status` |
| Stage | **+** | `git add <file>` / `git add .` |
| Commit | Message box → **Commit** | `git commit -m "msg"` |
| Push / pull | **Sync Changes** | `git push` / `git pull` |
| History | Timeline view, or *Git: View History* | `git log --oneline` |
| New branch | Status bar branch name | `git switch -c <name>` |
| Clone | *Git: Clone* | `git clone <url>` |

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| *Please tell me who you are* | Run the two `git config --global` lines in Part 0 |
| Push asks for a password, then fails | Passwords are not accepted: sign in via VS Code Accounts, or `gh auth login` |
| *rejected — fetch first* / *non-fast-forward* | Someone (or GitHub's web editor) changed the remote: `git pull`, resolve any conflict, then `git push` |
| Source Control panel is empty | You opened a single file, not the folder: *File → Open Folder* |
| *src refspec main does not match any* | No commits yet: commit before your first push |
