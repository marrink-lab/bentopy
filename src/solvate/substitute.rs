use eightyseven::structure::{Atom, AtomName};
use glam::Vec3;

use crate::{args::Substitute, placement::PlaceMap};

pub struct Substitution<'sol> {
    name: AtomName,
    placemap: PlaceMap<'sol>,
}

impl<'sol> Substitution<'sol> {
    pub fn new(name: impl Into<AtomName>, placemap: PlaceMap<'sol>) -> Self {
        Self {
            name: name.into(),
            placemap,
        }
    }

    pub fn natoms(&self) -> usize {
        self.placemap.unoccupied_count() as usize
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn iter_atoms(&self) -> impl Iterator<Item = Atom> + '_ {
        self.placemap.iter_positions().enumerate().map(|(i, pos)| {
            let num = i as u32 + 1;
            Atom {
                resnum: num,
                resname: self.name,
                atomname: self.name,
                atomnum: num,
                position: pos,
                velocity: Vec3::ZERO,
            }
        })
    }

    /// Fuse the contents of two [`Substitution`] maps.
    pub fn glue(&mut self, substitution: &Substitution<'sol>) {
        self.placemap &= substitution.placemap.clone();
    }
}

pub fn substitute<'sol>(
    rng: &mut rand::rngs::StdRng,
    placemap: &mut PlaceMap<'sol>,
    substitutes: &[Substitute],
) -> Vec<Substitution<'sol>> {
    use rand::seq::SliceRandom;

    // FIXME: Make this into a nicer error. Already explains the situation well, but the assert
    // failure is ugly in the output.
    let n_free = placemap.unoccupied_count();
    let n_subs = substitutes.iter().map(|s| s.number).sum::<u64>();
    assert!(
        n_subs <= n_free,
        "the number of substitutions that are specified ({n_subs}) exceeds the number of \
            solvent positions that can be substituted ({n_free})"
    );

    let mut substitutions = Vec::new();
    let mut zebra = Vec::new();
    for substitute in substitutes {
        eprint!(
            "\tSubstituting {} solvent locations with {}... ",
            substitute.number, substitute.name
        );
        let start = std::time::Instant::now();
        if substitute.number == 0 {
            eprintln!("Skipping.");
            continue;
        }

        let n_free = placemap.unoccupied_count();
        if n_free == 0 {
            eprintln!("No free solvent spots left.");
            break;
        }

        zebra.clear();
        zebra.resize(n_free as usize, false);
        zebra[..substitute.number as usize].fill(true);
        zebra.shuffle(rng);
        let mut zebra = zebra.iter();

        let mut old_placemap = placemap.clone(); // Becomes the map for our substitute.
        for mut placement in placemap.iter_placements_mut() {
            for i in 0..placement.len() {
                let is_occupied = placement.get(i);
                // If this position is not already occupied by a non-solvent bead, and if our
                // shuffled values indicate that we need to place a substitute here, we set this
                // bead to occupied in the solvent placemap.
                if !is_occupied && *zebra.next().unwrap() {
                    placement.set_occupied(i)
                }
            }
        }
        assert_eq!(zebra.len(), 0);

        old_placemap ^= placemap.clone();
        let substitute_placemap = !old_placemap;
        substitutions.push(Substitution::new(
            substitute.name.as_str(),
            substitute_placemap,
        ));
        eprintln!("Took {:.3} s.", start.elapsed().as_secs_f32());
    }

    substitutions
}
