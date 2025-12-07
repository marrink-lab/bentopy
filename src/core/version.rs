//! Provide additional version information, including the current git hash.

pub struct Version {
    pkg_version: &'static str,
    git_version: &'static str,
}

impl std::fmt::Display for Version {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self {
            pkg_version,
            git_version,
        } = self;
        write!(f, "{pkg_version} ({git_version})")
    }
}

// This makes it possible to use it as the version that Clap displays for a command.
impl Into<clap::builder::Str> for Version {
    fn into(self) -> clap::builder::Str {
        self.to_string().into()
    }
}

pub const VERSION: Version = Version {
    pkg_version: env!("CARGO_PKG_VERSION"),
    git_version: git_version::git_version!(
        args = ["--broken", "--always", "--exclude", "*"],
        prefix = "git:",
        fallback = "release"
    ),
};
