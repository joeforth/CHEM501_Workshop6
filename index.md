# Workshop 6 - Git, GitHub, and FAIR Practice

## Using Git from the Command Line

| --- | --- |
| Bring | Laptop, MKR, Nicla|
| Output | A project repo URL submitted on Canvas, containing code, README, and a history of commits |

## FAIR Practice


## Overview
In this Workshop, you'll create a GitHub repo that will host the code and documentation for your project. This is both a key part of making your work FAIR and 50% of your final grade. You'll start by *each* creating a GitHub repo and adding the basic components to get used to using git. You'll then work in your groups to create a GitHub repo for your project - each group should use and submit a single repo for their projects, and your repos will be assessed as a group.

Here we'll place special emphasis on getting used to running git in the command line / terminal in VS Code. This is both because it's useful and because John and Gabriella will rightly push you to use it in CHEM504 next semester. As with other Workshops, the hard bit is getting things set up - using the actual tools is the easy bit!

**How the workshop runs.** Every task is done with git commands typed into a terminal. The VS Code integrated terminal (opened with Ctrl+\` / Cmd+\`) is the most convenient, but any terminal works.

## Part 1 — Setup
Everyone needs git installed, a GitHub account, and a way for git to authenticate with GitHub before anything else works. You've probably done this already as part of Workshop 0, but if you haven't, the instructions are below. If something goes wrong, see the [Appendix](#appendix--setup-and-terminal-problems).

1. **Install git.** Windows: [git-scm.com](https://git-scm.com/downloads), accept the defaults (this installs Git Bash and Git Credential Manager). macOS: run `git --version` in a terminal and accept the Xcode tools prompt. If VS Code was open while you installed git, quit it completely and reopen it.
1. **Open VS Code and open a terminal.** Use *Terminal → New Terminal*.
    - **Windows:** Command Palette (Ctrl+Shift+P) → *Terminal: Select Default Profile* → **Git Bash**, then open a new terminal. This gives you the same commands as the Mac instructions; avoid *Command Prompt* (cmd).
    - **macOS:** Command Palette (Cmd+Shift+P) → *Shell Command: Install 'code' command in PATH*.
    - Run `git --version` in your terminal and check it prints a version number.
1. **Tell git who you are** (once per computer). Use the same email as you use with your GitHub account. The last line is useful but not essential - look up what it does.
    ```bash
    git config --global user.name "Joe Forth"
    git config --global user.email "j.forth@liverpool.ac.uk"
    git config --global init.defaultBranch main
    git config --global core.editor "code --wait"
    ```
1. **Sign VS Code in to GitHub.** Click the Accounts icon (person, bottom-left) → *Sign in with GitHub* → approve in the browser.

**Checkpoint:** `git --version` works in the VS Code terminal, `git config --global --list` shows your name and email, and the Accounts icon in VS Code shows your GitHub username.

## Part 2 — Five ideas you need
Git is a lab notebook for code: every commit is a dated, signed entry you can always return to; GitHub is the shared copy of that notebook online. Key bits of jargon:

1. **Repository (repo):** a project folder that git is tracking. The history lives in a hidden `.git` folder inside that folder. *Commands:* `git init` turns the current folder into a repo; `git clone <url>` copies an existing repo from the internet.
1. **Working → staging → committing:** you edit files (working), choose which changes go into the next entry (staging), then record them with a message (commiting). *Commands:* `git add <file>` stages a change; `git status` shows which files are in each state; `git diff` shows unstaged edits line by line, and `git diff --staged` shows what will go into the next commit.
1. **Commit:** Changes to the code that you're adding, plus a message saying *why*. Make small, frequent commits rather than one giant commit the night before the deadline. You're assessed on your commit record - both as evidence of incremental improvement and of teamwork. *Commands:* `git commit -m "message"` records everything staged; `git log --oneline` lists your commits.
1. **Remote:** the copy on GitHub, called `origin`. **Push** sends your commits up. **Pull** brings others' (or your other laptop's) commits down. *Commands:* `git clone <url>` sets `origin` up for you; `git push` and `git pull` do the sending and fetching.
1. **Branch:** a parallel line of work. `main` is the version that works - you can create new versions ('branches') to try risky changes that might not work. *Commands:* `git switch -c <name>` creates a branch and moves onto it; `git switch main` goes back; `git merge <name>` brings a finished branch into the one you're on.

The flow you should be repeating during your whole project: **edit → stage → commit → push**, or in commands: edit your files → `git add` → `git commit` → `git push`.

## Part 3 — Create your repo on GitHub
Before using git on your own computer, you'll create the repo on GitHub with the files every repo should start with.

In your browser, go to [github.com/new](https://github.com/new) and fill in:
- **Repository name:** e.g. `chem501-git-practice`.
- **Description:** one sentence saying what the repo is for. This is what makes it findable.
- **Visibility:** Private.
- **Add README:** on. The README is the front page of your repo, explaining what it is and how to use it.
- **Add .gitignore:** choose the **Python** template. This lists files git should never track (temporary and environment files); we'll add Arduino-specific lines in Part 5.
- **Choose a license:** **MIT License**. Without a licence, nobody else is legally allowed to reuse your code, however public it is.
- Click **Create repository**.

The description, README and licence are the starting point of FAIR practice: the description makes your work *findable*, the README tells others what it is and how to use it, and the licence makes it legally *reusable*. Every repo you make, including your group project repo, should have all three from the start.

**Checkpoint:** your new repo's page on GitHub shows `README.md`, `.gitignore` and `LICENSE`.

## Part 4 — Clone your repo and make your first commit (20 min)
Now you'll copy your repo onto your computer and make your first change with git from the command line.

**Task 4.1 — Clone your repo.** *Cloning* downloads a copy of a GitHub repo, with its full history, and connects it to GitHub so you can push and pull.

1. On your repo's GitHub page, click the green **Code** button, make sure **HTTPS** is selected, and copy the URL.
1. In the VS Code terminal, run `pwd` to see where you are. With no folder open, VS Code's terminal starts in your home folder, which is a good place to keep repos. (If VS Code has a folder open, use *File → Close Folder* first, and avoid folders synced by OneDrive or iCloud Drive, which clash with git.) Then run (pasting your own URL):
    ```bash
    git clone https://github.com/<your-username>/chem501-git-practice.git
    ```
    This creates a new folder called `chem501-git-practice` inside your current one. The first time, VS Code may ask you to sign in to GitHub; approve it.
1. Open the new folder: *File → Open Folder* → `chem501-git-practice`, then *Terminal → New Terminal*. Run `ls -a` and you should see `README.md`, `LICENSE`, `.gitignore` and the hidden `.git` folder.

**Task 4.2 — Make your first commit.**

1. Open `README.md` in VS Code. Under the title, add a sentence or two saying what the repo is for and who you are, and save it.
1. Check what git sees, then stage, commit and push:
    ```bash
    git status
    git add README.md
    git commit -m "Describe repo in README"
    git push
    ```
    Because you cloned the repo, git already knows where to push.

| What you did | Command |
| --- | --- |
| Copy a GitHub repo to your computer | `git clone <repo-url>` |
| See what's changed | `git status` |
| Stage a change | `git add <file>` |
| Commit | `git commit -m "message"` |
| Send it to GitHub | `git push` |

**Checkpoint:** Refresh your repo's page on GitHub; the README text should show your new sentence, and the commit count should have gone up by one.

## Part 5 — The everyday loop (25 min)

Most of your git use is one loop: edit, look at what changed, stage, commit, push.

As well as tracking what we do, we can also give git a list of files that we want it to ignore - this is important if you have files containing personal information (e.g., WiFi passwords).

**Adding to your `.gitignore`.** GitHub created a `.gitignore` for Python when you made the repo. Open it and add these lines at the top, so your WiFi credentials and editor clutter are never tracked:

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

**Task 5.1 — The loop.** Change the sensor read interval in your sketch, then:

```bash
git status                 # what has changed?
git diff                   # show the changes line by line
git add sensor_sketch.ino  # stage one file (or: git add .)
git commit -m "Change sensor read interval"
git push
git log --oneline          # see your history
```

**Task 5.2 — Again, from memory.** Add a comment block at the top of the sketch describing the board and sensor, and run through the same loop without looking at the commands above.

**Good commit messages** say what and why in under \~60 characters, written as a command:

| Weak | Better |
| --- | --- |
| `stuff` | `Add MQTT publish every 30 s` |
| `fixed it` | `Fix WiFi reconnect after router drop` |
| `final final v2` | `Calibrate temperature offset (−1.2 °C)` |

Aim for a commit each time something works or each lab session, whichever comes first.

## Part 6 — Branches, pulling, and undoing (25 min)

Branches let you try a change without breaking the version that works; pulling keeps two machines in step; a handful of commands undo most mistakes.

**Task 6.1 — Try a feature on a branch.** Add a second sensor reading (or change the upload format) on a new branch, merge it back once it works.

| Step | Command |
| --- | --- |
| New branch | `git switch -c add-humidity` |
| Commit and push work | `git add .` → `git commit -m "..."` → `git push -u origin add-humidity` |
| Merge into main | `git switch main` → `git merge add-humidity` → `git push` |

Alternatively, open a **pull request** on github.com and merge it there; this leaves a visible record of the change for your marker.

**Task 6.2 — Pull a change made elsewhere.** On github.com, edit `README.md` with the pencil icon and commit. Back in your terminal, run `git pull`. Always pull before you start work on a second computer.

**Task 6.3 — Undo.**

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

**Task 6.4:** Restructure your repo to the layout above, write a first README, commit and push. Check on github.com that it renders.

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
| *Please tell me who you are* | Run the two `git config --global` lines in Part 1 |
| Push asks for a password, then fails | Passwords are not accepted: sign VS Code in to GitHub and push from the VS Code terminal (see the Appendix) |
| *rejected — fetch first* / *non-fast-forward* | Someone (or GitHub's web editor) changed the remote: `git pull`, resolve any conflict, then `git push` |
| *fatal: not a git repository* | Your terminal is in the wrong folder: open the folder `git clone` created (*File → Open Folder*), or `cd` into it |

## Appendix — Setup and terminal problems

**Which shell on Windows?** Git commands work the same in Git Bash, PowerShell and Command Prompt (cmd), but the commands around them don't. In cmd, `ls` and `pwd` don't exist, and `#` isn't a comment, so pasting `git status   # what has changed?` passes the comment to git as arguments and fails. PowerShell mostly works. Git Bash matches the Mac instructions exactly, which is why we recommend it. The terminal's shell name appears in the terminal panel's top-right; the **+** dropdown next to it opens a different shell for a one-off.

**`git` is not recognized / command not found.** VS Code reads PATH only when it starts. Quit it completely (every window), reopen it, and open a new terminal. If it still fails on Windows, rerun the git installer and keep the default *Git from the command line and also from 3rd-party software*.

**Stuck in Vim.** If you see a screen full of `~` lines after `git commit`, `git merge` or `git pull`, git has opened the Vim editor. Press `Esc`, type `:wq`, press Enter to save and continue (or `:q!` to abort). Then run the `core.editor` line from Part 1 so it doesn't happen again. On a Mac, that line only works after *Install 'code' command in PATH*.

**`git push` asks for a username and password (common on Mac).** GitHub rejects account passwords. Sign VS Code in to GitHub (Part 1) and push from the VS Code integrated terminal, which uses that login as long as the setting `git.terminalAuthentication` is on (the default). Outside VS Code, install the [GitHub CLI](https://cli.github.com) and run `gh auth login`. On Windows, Git Credential Manager may instead open a browser window to approve.

**Repo inside OneDrive or iCloud.** Symptoms include files that are locked, duplicates with "conflicted copy" in the name, or errors about `.git/index.lock`. Move the project to an unsynced folder (`C:\Users\<you>\code`, `~/code`) and clone it fresh from GitHub.

**Cloned into the wrong place.** `git clone` creates the new folder inside whichever folder your terminal is in. If no folder is open in VS Code, that's your home folder. Run `pwd` before cloning. If the clone ends up somewhere you don't want, delete that folder and clone again: nothing is lost, because GitHub still has everything.

**Repo not found when cloning.** For a private repo, GitHub answers *Repository not found* when you're not signed in. Sign VS Code in to GitHub (Part 1) and clone from the VS Code terminal. Also check you copied the whole URL.

**Folder names with spaces.** Quote them: `cd "Workshop 6"`. Pressing Tab to autocomplete a name adds the quoting for you.

**Ctrl+\` doesn't open a terminal.** It doesn't work on some keyboard layouts. Use *Terminal → New Terminal* from the menu bar.

**`conda activate` fails in the VS Code terminal.** Each shell needs conda set up once: in Git Bash run `conda init bash`; in PowerShell run `conda init powershell` (if PowerShell then says *running scripts is disabled*, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`). Open a new terminal afterwards.

**Every line shows as changed in a mixed Mac/Windows team.** The cause is different line endings. Add a file called `.gitattributes` to the repo containing the single line `* text=auto`, then commit it.
