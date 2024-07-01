lf = LoddersFegley1998()
@test prod(length(x)==78 for x in lf )
@test lf.outer.Fe ≈ (lf.CI.Fe + lf.CK.Fe + lf.CM.Fe + lf.CO.Fe + lf.CR.Fe + lf.CV.Fe)/6
@test lf.inner.Fe ≈ (lf.H.Fe + lf.L.Fe + lf.LL.Fe + lf.R.Fe + lf.EH.Fe + lf.EL.Fe)/6
@test lf.inner.Rh ≈ (lf.H.Rh + lf.L.Rh)/2
@test lf.outer.Ta ≈ (lf.CI.Ta + lf.CM.Ta)/2