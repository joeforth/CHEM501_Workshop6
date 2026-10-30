# Workshop 6 — Git, GitHub, and FAIR Practice

## FAIR Practice
FAIR is the act of making your code and data Findable, Accessible, Interoperable, and Reusable. As a formal idea-with-a-name, it's been around for about a decade, with the main paper in the area being [Wilkinson et al.](https://www.nature.com/articles/sdata201618). The idea is that we should be making our code and data as FAIR as possible, while also respecting things like commercial sensitivity. Today, we'll start the process of making the code in your CHEM501 mini-projects FAIR by creating a GitHub repo (making it findable and accessible), writing a README using Markdown (documenting it to make it interoperable, giving it a licence to make it reusable), and managing the code in a repo with git. 

You'll start by *each* creating a GitHub repo and adding the basic components to get used to using git. You'll then work in your groups to create a GitHub repo for your project - each group should use and submit a single repo for their projects, and your repos will be assessed as a group. Your GitHub repo constitutes 50% of the grade for CHEM501.

Here we'll place special emphasis on getting used to running git in the command line / terminal in VS Code. This is both because it's useful and because you will (rightly) be pushed to use it this way in CHEM504 next semester. As with other Workshops, the hard bit is getting things set up - using the actual tools is the easy bit!

**How the workshop runs.** Apart from a few steps on the GitHub website, every task is done with git commands typed into a terminal. The VS Code integrated terminal (opened with Ctrl+\` / Cmd+\` or by clicking Terminal → New Terminal) is preferable, but any terminal works.

## Five Ideas You Need
Git is a lab notebook for code: every commit is a dated entry, with your name on it, that you can always return to; GitHub is the shared copy of that notebook online. Key bits of jargon:

1. **Repository (repo):** a project folder that git is tracking. The history lives in a hidden `.git` folder inside that folder. *Commands:* `git init` turns the current folder into a repo; `git clone <url>` copies an existing repo from the internet.
1. **Working → staging → committing:** you edit files (working), choose which changes go into the next entry (staging), then record them with a message (committing). *Commands:* `git add <file>` stages a change; `git status` shows which files are in each state; `git diff` shows unstaged edits line by line, and `git diff --staged` shows what will go into the next commit.
1. **Commit:** a set of changes to the code that you're adding, plus a message saying *why*. Make small, frequent commits rather than one giant commit the night before the deadline. You're assessed on your commit record - both as evidence of incremental improvement and of teamwork. *Commands:* `git commit -m "message"` records everything staged; `git log --oneline` lists your commits.
1. **Remote:** the copy on GitHub, called `origin`. **Push** sends your commits up. **Pull** brings others' (or your other laptop's) commits down. *Commands:* `git clone <url>` sets `origin` up for you; `git push` and `git pull` do the sending and fetching.
1. **Branch:** a parallel line of work. `main` is the version that works - you can create new versions ('branches') to try risky changes that might not work. *Commands:* `git switch -c <name>` creates a branch and moves onto it; `git switch main` goes back; `git merge <name>` brings a finished branch into the one you're on.

The flow you should be repeating during your whole project: **edit → stage → commit → push**, or in commands: edit your files → `git add` → `git commit` → `git push`.

## Part 1 — Setup
You need to set up git, GitHub, and VS Code so it works with git. You've likely done at least some of this as part of the setup before CHEM501 started, but some of the setup for this workshop will be new. Work through the steps below to get things set up. If something goes wrong, see the [Appendix](#appendix--setup-and-terminal-problems).

1. **Install git.** Windows: [git-scm.com](https://git-scm.com/downloads), accept the defaults (this installs Git Bash and Git Credential Manager). macOS: run `git --version` in a terminal and accept the Xcode tools prompt. If VS Code was open while you installed git, quit it completely and reopen it.
1. **Open VS Code and open a terminal.** Use *Terminal → New Terminal*.
    - **Windows:** Command Palette (Ctrl+Shift+P) → *Terminal: Select Default Profile* → **Git Bash**, then open a new terminal. This gives you the same commands as the Mac instructions; avoid *Command Prompt* (cmd).
    - **macOS:** Command Palette (Cmd+Shift+P) → *Shell Command: Install 'code' command in PATH*.
    - Run `git --version` in your terminal and check it prints a version number.
1. **Tell git who you are** (once per computer). Use your own name, and the same email as you use with your GitHub account. The last two lines are useful but not essential - look up what they do.
    ```bash
    git config --global user.name "Your Name"
    git config --global user.email "you@liverpool.ac.uk"
    git config --global init.defaultBranch main
    git config --global core.editor "code --wait"
    ```
1. **Sign VS Code in to GitHub.** Click the Accounts icon (person, bottom-left) → *Sign in with GitHub* → approve in the browser.

**Checkpoint:** `git --version` works in the VS Code terminal, `git config --global --list` shows your name and email, and the Accounts icon in VS Code shows your GitHub username.

## Part 2 — Create your repo on GitHub
Before using git on your own computer, you'll create the repo on GitHub with the files every repo should start with.

In your browser, go to [github.com/new](https://github.com/new) and fill in:
- **Repository name:** e.g. `chem501-git-practice`.
- **Description:** one sentence saying what the repo is for.
- **Visibility:** Public (otherwise it's not Findable!). Your group project repo must be public too.
- **Add README:** On. The README is the front page of your repo, explaining what it is and how to use it.
- **Add .gitignore:** choose the **Python** template. This lists files git should never track (temporary and environment files).
- **Choose a license:** Up to you - the MIT License is a good catch-all option, but most of them are suitable.
- Click **Create repository**.

The description, README, and licence are the starting point of FAIR practice: the description makes your work *findable*, the README tells others what it is and how to use it, and the licence makes it legally *reusable*. Every repo you make, including your group project repo, should have all three from the start.

**Checkpoint:** Check your new repo's page on GitHub shows `README.md`, `.gitignore`, and `LICENSE`.

## Part 3 — Clone your repo and make your first commit
Now you'll copy your repo onto your computer and make your first change with git from the command line.

**Task 3.1 — Clone your repo.** *Cloning* downloads a copy of a GitHub repo, with its full history, and connects it to GitHub so you can push and pull.

1. On your repo's GitHub page, click the green **Code** button, make sure **HTTPS** is selected, and copy the URL.
1. In the VS Code terminal, run `pwd` to see where you are. With no folder open, VS Code's terminal starts in your home folder, which is a good place to keep repos. (If VS Code has a folder open, use *File → Close Folder* first, and avoid folders synced by OneDrive or iCloud Drive, which clash with git.) Then run (pasting your own URL):
    ```bash
    git clone https://github.com/<your-username>/chem501-git-practice.git
    ```
    This creates a new folder called `chem501-git-practice` inside your current one. The first time, VS Code may ask you to sign in to GitHub; approve it.
1. Open the new folder: *File → Open Folder* → `chem501-git-practice`, then *Terminal → New Terminal*. Run `ls -a` and you should see `README.md`, `LICENSE`, `.gitignore` and the hidden `.git` folder.

**Task 3.2 — Make your first commit.**

1. Open `README.md` in VS Code. README.md docs are written in a language called [Markdown](https://www.markdownguide.org/). Use the [Markdown cheat sheet](https://www.markdownguide.org/cheat-sheet/) to add a heading, sub-heading, and a numbered list to your README.md. Save the edits.
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

**Checkpoint:** Refresh your repo's page on GitHub; the README should show your new heading, sub-heading and list, and the commit count should have gone up by one.

## Part 4 — The everyday loop

Most of your git use is one loop: edit, look at what changed, stage, commit, push. We'll practise it here.

**Task 4.1 — Keep your passwords out of git.** Anything you push to GitHub stays in the history, even if you delete it later, so passwords must never be committed. The Arduino sketches therefore keep the WiFi name and password in a separate file, `arduino_secrets.h`, in the same folder as the `.ino` file (the Arduino IDE shows it as its own tab). It looks like this:

```c
#define SECRET_SSID "your-network-name"
#define SECRET_PASS "your-network-password"
```

The sketch loads it with `#include "arduino_secrets.h"`. Because the secrets are in their own file, you can tell git never to track that file, using `.gitignore`. GitHub created a `.gitignore` for Python when you made the repo. Open it and add these lines at the top:

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

Commit and push `.gitignore` *before* you add any Arduino code: `git add .gitignore` → `git commit -m "Ignore WiFi secrets and editor files"` → `git push`. Anyone who clones your repo will need to create their own `arduino_secrets.h`, so say how in your README.

**Task 4.2 — The loop.** Copy the Arduino and Python code that you used earlier in the workshop into your repo folder, then run:

```bash
git status                 # what has changed? arduino_secrets.h should NOT be listed
git diff                   # show changes to tracked files line by line
git add .                  # stage everything (or: git add <file> for one file)
git commit -m "Add MQTT sketch and Python notebook"
git push
git log --oneline          # see your history
```

Always give each commit a message with `-m`. If you forget, git opens an editor for you to write one: VS Code, if you ran the `core.editor` line in Part 1 (type the message, save, close the tab). If you see a screen full of `~` lines instead, you're in Vim; see the [Appendix](#appendix--setup-and-terminal-problems).

**Good commit messages** say what and why in under \~60 characters, written as a command:

| Bad Commit Message | Better Commit Message |
| --- | --- |
| `stuff` | `Add MQTT publish every 30 s` |
| `fixed it` | `Fix WiFi reconnect after router drop` |
| `final final v2` | `Calibrate temperature offset (−1.2 °C)` |

Aim for a commit each time something works or each lab session, whichever comes first.

## Part 5 — Branches, pulling, and undoing

Branches let you try making a change to your code without breaking the version that works; pulling keeps two machines in step; a handful of commands undo most mistakes. As you'll be working collaboratively on your code, it's essential you learn how to merge and pull code.

**Task 5.1 — Try a feature on a branch.** Modify the Arduino code, so that it publishes messages more frequently.

| Step | Command |
| --- | --- |
| New branch | `git switch -c faster-publish` |
| Commit and push work | `git add .` → `git commit -m "..."` → `git push -u origin faster-publish` |
| Merge into main | `git switch main` → `git merge faster-publish` → `git push` |

Alternatively, merge using a **pull request** (PR) on github.com. A pull request is a request to merge one branch into another. GitHub shows it as a page with the changes line by line, so you (or a teammate) can review them, comment on them, and push more commits to fix problems before anything reaches `main`. In team projects this is the usual way to merge: someone else checks your work before it goes in.

To open a pull request, push your branch as above, then on your repository's page on github.com click **Compare & pull request** (or go to the **Pull requests** tab → **New pull request**, and choose `faster-publish` to merge into `main`). Add a title and a short description of what you changed and why, then click **Create pull request**. When you're happy with it, click **Merge pull request**. Finally, update your local copy with `git switch main` → `git pull`.

Your commits appear in the history whichever way you merge, but a pull request also records the branch, the changes and the merge in one place. That makes it easy for your marker to see that you used a branch.

**Task 5.2 — Pull a change made elsewhere.** On github.com, edit `README.md` with the pencil icon and commit. Back in your terminal, run `git pull`. Always pull before you start work on a second computer.

**Task 5.3 — Undo.**

| Situation | Command |
| --- | --- |
| Throw away uncommitted edits to a file | `git restore <file>` |
| Unstage a file | `git restore --staged <file>` |
| Fix the last commit message (not yet pushed) | `git commit --amend -m "Better message"` |
| Undo a pushed commit safely | Find the commit with `git log --oneline`, then `git revert <commit-id>` |

**Merge conflicts.** If you and GitHub both changed the same lines, `git pull` stops and reports a conflict. `git status` lists the affected files; inside each, git marks the clash between `<<<<<<<`, `=======` and `>>>>>>>` lines. Edit the file to keep the version you want, delete the markers, then `git add` and `git commit`. Pulling before you start each session avoids most conflicts.

## Part 6 — Start Your Projects
In your pairs / groups / whatever - start building the GitHub repo for your projects. You're assessed in groups, so you should all work on a single repo, which means you all need access to it.

1. One of you creates the project repo the same way as your practice repo (Part 2): public, with a README, the Python `.gitignore` and a licence. Clone it and add the lines from Task 4.1 to `.gitignore` straight away.
2. On the repo's page on github.com, go to **Settings → Collaborators** (GitHub may ask you to confirm your password).
3. Click **Add people**, search for each group member's GitHub username, and click **Add to repository**. Add your module lead the same way.
4. Each person then accepts the invitation, either from the email GitHub sends them or at github.com/notifications. The invite expires after 7 days.

Once they've accepted, collaborators can push directly to the repo. Everyone else clones it with `git clone <url>`. Agree as a group whether you'll push straight to `main` or work on branches and merge with pull requests (Part 5).

**README checklist** (write it in Markdown; VS Code previews with Ctrl+Shift+V):

- [ ] One-paragraph summary: what it measures, and why
- [ ] Hardware list and a wiring diagram or photo
- [ ] Libraries needed, with versions
- [ ] How to build and upload, including creating `arduino_secrets.h`
- [ ] Where the data goes (serial, cloud dashboard, MQTT) and a sample of it
- [ ] Known limitations and what you would do next

Now - start your projects.

## Before you submit your final project repo

- [ ] Repo is public (so it's Findable and Accessible) and your module lead's GitHub account is added under *Settings → Collaborators*
- [ ] `git status` says *nothing to commit, working tree clean* and everything is pushed
- [ ] Code compiles from a fresh clone after creating `arduino_secrets.h`
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
| Change branch | `git switch <name>` |
| Merge a branch into the current one | `git merge <name>` |
| Throw away uncommitted edits | `git restore <file>` |
| Clone | `git clone <url>` |

## Troubleshooting

| Symptom | Fix |
| --- | --- |
| *Please tell me who you are* | Run the `user.name` and `user.email` lines from Part 1 |
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

**Repo not found when cloning.** Usually the URL is wrong or incomplete: copy it again from the green **Code** button. For a private repo (someone else's, for example), GitHub also answers *Repository not found* when you're not signed in. Sign VS Code in to GitHub (Part 1) and clone from the VS Code terminal.

**Folder names with spaces.** Quote them: `cd "Workshop 6"`. Pressing Tab to autocomplete a name adds the quoting for you.

**Ctrl+\` doesn't open a terminal.** It doesn't work on some keyboard layouts. Use *Terminal → New Terminal* from the menu bar.

**`conda activate` fails in the VS Code terminal.** Each shell needs conda set up once: in Git Bash run `conda init bash`; in PowerShell run `conda init powershell` (if PowerShell then says *running scripts is disabled*, run `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`). Open a new terminal afterwards.

**Every line shows as changed in a mixed Mac/Windows team.** The cause is different line endings. Add a file called `.gitattributes` to the repo containing the single line `* text=auto`, then commit it.
