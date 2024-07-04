use std::str::FromStr;

use glam::Vec3;

#[derive(Debug, Default, Clone)]
pub struct Limits {
    pub minx: Option<f32>,
    pub maxx: Option<f32>,
    pub miny: Option<f32>,
    pub maxy: Option<f32>,
    pub minz: Option<f32>,
    pub maxz: Option<f32>,
}

impl Limits {
    pub fn is_inside(&self, point: impl Into<Vec3>) -> bool {
        let point = point.into();
        self.minx.map_or(true, |minx| minx <= point.x)
            && self.miny.map_or(true, |miny| miny <= point.y)
            && self.minz.map_or(true, |minz| minz <= point.z)
            && self.maxx.map_or(true, |maxx| point.x <= maxx)
            && self.maxy.map_or(true, |maxy| point.y <= maxy)
            && self.maxz.map_or(true, |maxz| point.z <= maxz)
    }

    // NOTE: This implementation is slightly verbose, and could be shortened. But that is not
    // necessary. This is not speed-critical and I'd rather keep it entirely clear what is
    // going on.
    /// Determine the appropriate box size based on the constraints posed by the [`Limits`].
    pub fn box_size(&self, boxsize: Vec3) -> Vec3 {
        // We bake the limits with the information from the original box size.
        let baked_limits = Limits {
            minx: Some(self.minx.unwrap_or(0.0)),
            maxx: Some(self.maxx.unwrap_or(boxsize.x)),
            miny: Some(self.miny.unwrap_or(0.0)),
            maxy: Some(self.maxy.unwrap_or(boxsize.y)),
            minz: Some(self.minz.unwrap_or(0.0)),
            maxz: Some(self.maxz.unwrap_or(boxsize.z)),
        };

        // We can safely unwrap these, because we just initialized them as Some(f32).
        Vec3::new(
            baked_limits.maxx.unwrap() - baked_limits.minx.unwrap(),
            baked_limits.maxy.unwrap() - baked_limits.miny.unwrap(),
            baked_limits.maxz.unwrap() - baked_limits.minz.unwrap(),
        )
    }
}

impl FromStr for Limits {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let Ok([minx, maxx, miny, maxy, minz, maxz]): Result<[Option<f32>; 6], _> = s
            .split(',')
            .map(|v| v.trim().parse().ok())
            .collect::<Vec<_>>()
            .try_into()
        else {
            let n = s.split(',').count();
            return Err(format!("expected 6 values, found {n}"));
        };

        Ok(Self {
            minx,
            maxx,
            miny,
            maxy,
            minz,
            maxz,
        })
    }
}
