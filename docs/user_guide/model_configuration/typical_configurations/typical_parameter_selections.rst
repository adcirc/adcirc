.. _typical_adcirc_parameter_selections:

Typical Parameter Selections
============================

.. meta::
   :description: Typical Adcirc Parameter Selections in ADCIRC
   :keywords: adcirc, typical, adcirc, parameter, selections

**Parameter selections** concern the model parameter choices one must make when setting up and running simulations. ADCIRC and the coupled ADCRC+SWAN model have a large number of user-controllable parameters. Some parameters represent a balance in model accuracy vs. stability, like enabling vs. disabling advection via :ref:`NOLICA` and :ref:`NOLICAT`. Others balance accuracy vs. run time, like the implicit vs. lumped explicit GWCE formulation via :ref:`IM` (this also affects stability). Meanwhile some parameters represent values for which the true "best" choice is case-dependent or the subject of ongoing research, like controls affecting drag formulations at the sea surface and sea floor. This page is meant to provide users with example model formulations used by experienced ADCIRCers, and to provide a venue for discussion of parameter selection.

In Storm Surge
--------------

The following table presents a variety of parameters as used in storm surge modeling by various groups. Parameters in this table cover a range of model setup choices, and include parameterizations built for forecasting-focused simulations, storm hazard-focused simulations, and hindcast-focused simulations. Note that in all cases, nodal attributes `Manning's n at sea floor </Manning%27s_n_at_sea_floor "Manning's n at sea floor">`_ and `surface directional effective roughness length </Surface_directional_effective_roughness_length "Surface directional effective roughness length">`_ are in use.

.. list-table:: Storm Surge Parameter Selections
   :widths: 10 15 8 10 8 10 10 8 8 8 5 8 8 5 5 15 8 8 8 8
   :header-rows: 1
   :stub-columns: 1
   :class: sticky-table

   * - Region
     - Mesh
     - np
     - GWCE Formulation
     - Wind Drag Law
     - Upper Wind Drag Limit
     - Min Bottom Drag Coef
     - Advective Terms
     - Steric
     - Vert. Datum
     - ESL
     - ESL Node Count
     - time step
     - h0
     - velmin
     - tidal constituents
     - convcr
     - SWAN time step
     - SWAN MXITNS
     - SWAN NPNTS
   * - Louisiana
     - LA_v17a-WithUpperAtch_chk.grd
     - 1593521
     - implicit
     - Powell
     - 0.002
     - 0.0
     - off
     - 0.228184
     - navd
     - 0.05
     - all
     - 1.0
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - Texas
     - tx2008_r35h.grd
     - 3352598
     - implicit
     - Garratt
     - 0.002
     - 0.0
     - off
     - 0.276300
     - navd?
     - none
     - none
     - 1.0
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - National
     - hsofs.14
     - 1813443
     - implicit or explicit
     - Garratt
     - 0.0028
     - 0.0025
     - off
     - none
     - msl
     - none
     - none
     - 2.0
     - 0.05
     - 0.05
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - National (Drag Experiment)
     - hsofs.14
     - 1813443
     - explicit
     - Garratt
     - 0.002
     - 0.0
     - off
     - none
     - msl
     - none
     - none
     - 2.0
     - 0.05
     - 0.05
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - MS/AL and FL Panhandle
     - NGOM_RT_v18j_chk.grd
     - 2051346
     - implicit
     - Powell
     - 0.002
     - 0.0
     - off
     - 0.230000
     - navd
     - 0.02
     - all
     - 1.0
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - North Carolina (low res)
     - nc_inundation_v6d_rivers_msl.grd
     - 295328
     - explicit
     - Garratt
     - 0.0035
     - 0.003
     - on
     - none
     - msl
     - none
     - none
     - 0.5
     - 0.02
     - 0.02
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-10
     - 1200
     - 20
     - 95
   * - North Carolina (high res)
     - nc_inundation_v9.99a_w_rivers.grd
     - 624782
     - explicit
     - Garratt
     - 0.0028
     - 0.003
     - on
     - 0.0
     - msl
     - none
     - none
     - 0.5
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-10
     - 1200
     - 20
     - 95
   * - Delmarva
     - FEMA_R3_20110303_MSL.grd
     - 1875689
     - explicit
     - Garratt
     - 0.0035
     - 0.003
     - on
     - none
     - msl
     - none
     - none
     - 1.0
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - 1200
     - 20
     - 95
   * - NY/NJ
     - FEMA_R2_norivers_gcs_mNAVD.grd
     - 604790
     - implicit
     - Garratt
     - 0.0035
     - 0.003
     - on
     - none
     - msl
     - none
     - none
     - 1.0
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,q1
     - 1.00E-08
     - 1200
     - 10
     - 95
   * - New England
     - NAC2014_R01_ClosedRivers.grd
     - 3110470
     - implicit
     - Powell
     - 0.002
     - 0.0
     - off
     - 0.10900
     - msl
     - none
     - none
     - ?
     - 0.1
     - 0.01
     - m2,s2,n2,k1,k2,o1,p1,q1
     - 1.00E-07
     - n/a
     - n/a
     - n/a

In Tides
--------

Many tide models are configured the same as the above surge models, although some users have reported poor tide performance if the minimum bottom drag coefficient is set to zero. One explanation of this difference in behavior is that the bottom drag differs between tidal and wind-driven flows. The vertical variation in horizontal flow is not known in a 2D model, making prescription of bottom drag ambiguous because the near-bottom velocity is unknown.

Discouraged Parameter Selections
--------------------------------

As ADCIRC has grown, some features have been improved upon, making others obsolete. This section addresses such parameters and values.

`NWS=19` should not be used, as `NWS=20` has the same input requirements but is considered to produce a better, more physically representative wind field. See :ref:`nws_parameter` and :ref:`gahm` for more details.
