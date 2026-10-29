# Workshop 6 - Git, GitHub, and FAIR Practice
## Using Git from the Command Line

|  |  |
| --- | --- 
| Bring | Laptop, MKR, Nicla|
| Output | A project repo URL submitted on Canvas, containing code, README, and a history of commits |

## FAIR Practice


## Overview
In this Workshop, you'll create a GitHub repo for your projects. You'll start by *each* creating a GitHub repo and adding the basic components to get used to using Git. You'll then work in your groups to create a GitHub repo for your project - each group should use a single repo for their projects, and your repos will be assessed as a group.

Special emphasis is place on getting used to running Git in the command line / terminal in VS Code. As with other Workshops, the hard bit is getting things setup - using the actual tools is the easy bit!

**How the workshop runs.** Every task is done with git commands typed into a terminal. The VS Code integrated terminal (opened with Ctrl+\` / Cmd+\`) is the most convenient, but any terminal works.

## Part 0 — Setup (15 min)

Everyone needs git installed, a GitHub account, and a way for git to authenticate with GitHub before anything else works. You've hopefully done this already as part of Workshop 0, but if you haven't, the instructions are below. If something goes wrong, see the [Appendix](#appendix--setup-and-terminal-problems).

1. **Install git.** Windows: [git-scm.com](https://git-scm.com/downloads), accept the defaults (this installs Git Bash and Git Credential Manager). macOS: run `git --version` in a terminal and accept the Xcode tools prompt. If VS Code was open while you installed git, quit it completely and reopen it.
1. **Open VS Code and open a terminal.** Use *File → Open Folder* to open your project folder, then *Terminal → New Terminal*.
    - **Windows:** Command Palette (Ctrl+Shift+P) → *Terminal: Select Default Profile* → **Git Bash**, then open a new terminal. This gives you the same commands as the Mac instructions; avoid *Command Prompt* (cmd).
    - **macOS:** Command Palette (Cmd+Shift+P) → *Shell Command: Install 'code' command in PATH*.
    - Check that `git --version` prints a version number.
1. **Tell git who you are** (once per computer). Use the same email as your GitHub account. The last line makes git open VS Code, rather than Vim, when it needs you to write a message:
```bash
git config --global user.name "Ada Lovelace"
git config --global user.email "ada@liverpool.ac.uk"
git config --global init.defaultBranch main
git config --global core.editor "code --wait"
```
1. **Sign VS Code in to GitHub.** Click the Accounts icon (person, bottom-left) → *Sign in with GitHub* → approve in the browser. GitHub doesn't accept passwords over the command line; once VS Code is signed in, `git push` from the VS Code terminal uses that login instead.
1. **Keep projects out of synced folders.** OneDrive (often your Windows Desktop and Documents) and iCloud Drive clash with git. Use a folder such as `C:\Users\<you>\code` or `~/code`.

**Checkpoint:** `git --version` works in the VS Code terminal, `git config --global --list` shows your name and email, and the Accounts icon shows your GitHub username.

## Part 1 — Five ideas you need (10 min)
Git is a lab notebook for code: every commit is a dated, signed entry you can always return to; GitHub is the shared copy of that notebook online.

1. **Repository (repo):** a project folder that git is tracking. The history lives in a hidden `.git` folder inside it. *Commands:* `git init` turns the current folder into a repo; `git clone <url>` copies an existing one from GitHub.
2. **Working → staged → committed:** you edit files (working), choose which changes go into the next entry (staged), then record them with a message (commit). *Commands:* `git add <file>` stages a change; `git status` shows which files are in each state; `git diff` shows unstaged edits line by line, and `git diff --staged` shows what will go into the next commit.
3. **Commit:** Changes to the code that you're adding, plus a message saying *why*. Make small, frequent commits rather than one giant commit the night before the deadline. You're assessed on your commit record - both as evidence of incremental improvement and of teamwork. *Commands:* `git commit -m "message"` records everything staged; `git log --oneline` lists your commits.
4. **Remote:** the copy on GitHub, called `origin`. **Push** sends your commits up. **Pull** brings others' (or your other laptop's) commits down. *Commands:* `git remote add origin <url>` connects a repo to GitHub (once); `git push` and `git pull` do the sending and fetching.
5. **Branch:** a parallel line of work. `main` is the version that works; try risky changes on a separate branch. *Commands:* `git switch -c <name>` creates a branch and moves onto it; `git switch main` goes back; `git merge <name>` brings a finished branch into the one you're on.

The flow you should be repeating during your whole project: **edit → stage → commit → push**, or in commands: edit your files → `git add` → `git commit` → `git push`.

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
| Push asks for a password, then fails | Passwords are not accepted: sign VS Code in to GitHub and push from the VS Code terminal (see the Appendix) |
| *rejected — fetch first* / *non-fast-forward* | Someone (or GitHub's web editor) changed the remote: `git pull`, resolve any conflict, then `git push` |
| *fatal: not a git repository* | Your terminal is in the wrong folder: `cd` into your project folder (or run `git init` if you haven't yet) |
| *src refspec main does not match any* | No commits yet: commit before your first push |

## Appendix — Setup and terminal problems

**Which shell on Windows?** Git commands work the same in Git Bash, PowerShell and Command Prompt (cmd), but the commands around them don't. In cmd, `ls` and `pwd` don't exist, and `#` isn't a comment, so pasting `git status   # what has changed?` passes the comment to git as arguments and fails. PowerShell mostly works. Git Bash matches the Mac instructions exactly, which is why we recommend it. The terminal's shell name appears in the terminal panel's top-right; the **+** dropdown next to it opens a different shell for a one-off.

**`git` is not recognized / command not found.** VS Code reads PATH only when it starts. Quit it completely (every window), reopen it, and open a new terminal. If it still fails on Windows, rerun the git installer and keep the default *Git from the command line and also from 3rd-party software*.

**Stuck in Vim.** If you see a screen full of `~` lines after `git commit`, `git merge` or `git pull`, git has opened the Vim editor. Press `Esc`, type `:wq`, press Enter to save and continue (or `:q!` to abort). Then run the `core.editor` line from Part 0 so it doesn't happen again. On a Mac, that line only works after *Install 'code' command in PATH*.

**`git push` asks for a username and password (common on Mac).** GitHub rejects account passwords. Sign VS Code in to GitHub (Part 0) and push from the VS Code integrated terminal, which uses that login as long as the setting `git.terminalAuthentication` is on (the default). Outside VS Code, install the [GitHub CLI](https://cli.github.com) and run `gh auth login`. On Windows, Git Credential Manager may instead open a browser window to approve.

**Repo inside OneDrive or iCloud.** Symptoms include files that are locked, duplicates with "conflicted copy" in the name, or errors about `.git/index.lock`. Move the project to an unsynced folder (`C:\Users\<you>\code`, `~/code`) and clone it fresh from GitHub.

**`git init` in the wrong place.** If no folder is open in VS Code, the terminal starts in your home folder, and `git init` there makes your whole home folder a repo. Run `pwd` before `git init`. To undo a mistaken `git init`, delete the hidden `.git` folder it created (and only that folder).

**Folder names with spaces.** Quote them: `cd "Workshop 6"`. Pressing Tab to autocomplete a name adds the quoting for you.

**Ctrl+\` doesn't open a terminal.** It doesn't work on some keyboard layouts. Use *Terminal → New Terminal* from the menu bar.

**`conda activate` fails in the VS Code terminal.** Each shell needs conda set up once: in Git Bash run `conda init bash`; in PowerShell run `conda init powershell` (if PowerShell then says *running scripts is disabled*, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`). Open a new terminal afterwards.

**Every line shows as changed in a mixed Mac/Windows team.** The cause is different line endings. Add a file called `.gitattributes` to the repo containing the single line `* text=auto`, then commit it.
