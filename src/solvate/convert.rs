/// Convert between `Vec3` types.
///
/// In this program we use both [`glam::Vec3`] and the [`eightyseven::structure::Vec3`] type, which
/// is an export of some version of `glam::Vec3` used in `eightyseven`.
///
/// This is very annoying, because we can't just use them interchangeably despite the types being
/// otherwise identical.
pub trait Convert<V> {
    /// Convert between `Vec3` types.
    fn convert(&self) -> V;
}

impl Convert<glam::Vec3> for eightyseven::structure::Vec3 {
    #[inline(always)]
    fn convert(&self) -> glam::Vec3 {
        glam::Vec3::from_array(self.to_array())
    }
}

impl Convert<eightyseven::structure::Vec3> for glam::Vec3 {
    #[inline(always)]
    fn convert(&self) -> eightyseven::structure::Vec3 {
        eightyseven::structure::Vec3::from_array(self.to_array())
    }
}
