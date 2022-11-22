creds <- git2r::cred_env("DOCKER_GITLAB_USER", "DOCKER_GITLAB_PW")
devtools::install_git("https://git.ecdf.ed.ac.uk/oalmelid/BaalChIP.git", credentials = creds, ref = "dev")
