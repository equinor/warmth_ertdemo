import numpy as np
from temperer.forward_modelling import Forward_model

def test_pad_multirift_results():
    pad_func = Forward_model._equalise_array_shape
    holder =  np.zeros((10, 2),dtype=float)
    new =  np.zeros((15, 2),dtype=float)
    r,_ = pad_func(holder,new,np.nan)
    assert r.shape == new.shape
    assert np.all(np.isnan(r[10:,:]))
    holder =  np.zeros((15, 2),dtype=float)
    new =  np.zeros((10, 2),dtype=float)
    _,r = pad_func(holder,new,np.nan)
    assert r.shape == holder.shape
    assert np.all(np.isnan(r[10:,:]))

# def instantiate_test_model():
#     example_model_path = './tests/data/example_model.pickle'
#     with open(example_model_path, 'rb') as file:
#         model = pickle.load(file)
#     f_model = forward_modelling.Forward_model(
#         model.parameters, model.builder.all_nodes[0])
#     return f_model


# def test_sediment_density():
#     f_model = instantiate_test_model()
#     f_model._parameters.rhowater = 1000
#     result = f_model._sediment_density(
#         np.array([0.2, 0.5]), np.array([2600, 2700]))
#     assert np.allclose(result, np.array([2280., 1850.])) == True


# def test_sediment_conductivity():
#     result = forward_modelling.Forward_model._sediment_conductivity(
#         np.array([0.2, 0.5]), np.array([2, 3]))
#     assert np.allclose(result, np.array([1.57200617, 1.34164079])) == True


# def test_check_beta():
#     f_model = instantiate_test_model()
#     results = f_model._check_beta(
#         100, 2.5, np.array([1.1, 2]), np.array([-2000, -1000]))
#     assert results[0] == True
#     results = f_model._check_beta(-100, 2.5,
#                                   np.array([1.1, 2]), np.array([-2000, -1000]))
#     assert results[0] == False


# def test_check_convergence():
#     f_model = instantiate_test_model()
#     T_this_step = np.array([10, 100, 300, 500])
#     T_last_step = np.array([8, 80, 150, 300])
#     f_model._parameters.convergence = 5
#     result = f_model._check_convergence(T_last_step, T_this_step)
#     assert result == False
#     f_model._parameters.convergence = 6
#     result = f_model._check_convergence(T_last_step, T_this_step)
#     assert result == True


# def test_build_crust_lithosphere_properties():
#     coord = np.array([0, 100, 200, 300, 400, 500, 600])
#     result = forward_modelling.Forward_model._build_crust_lithosphere_properties(
#         coord, 200, 400, 1000, 2000, 3000)
#     assert np.allclose(result, np.array(
#         [1000., 1000., 2000., 2000., 3000., 3000.])) == True


# def test_effective_density():
#     f_model = instantiate_test_model()
#     T_arr = np.array([10, 50, 100])
#     density = np.array([2000, 2500])
#     f_model._parameters.alphav = 3.0e-5
#     f_model.current_node.T0 = 5
#     result = f_model._effective_density(density, T_arr)
#     assert np.allclose(result, np.array([[1998.5, 2494.75]])) == True


# def test_generate_mesh_1D():
#     f_model = instantiate_test_model()
#     f_model._parameters.resolution = 100
#     f_model.current_node.hc = 150
#     f_model.current_node.hLith = 450
#     f_model._generate_mesh_1D()
#     assert np.allclose(f_model.current_node.coord_initial[:12], np.array(
#         [0.0, 100.0, 150.0, 200.0, 300.0, 400.0, 450.0, 500.0, 600.0, 700.0, 800.0, 900.0])) == True
#     assert f_model.current_node.ncrust == 3


# def test_remesh_crust_lith_asth():
#     coord = np.array([0, 100, 200, 300, 400, 500, 600])
#     result = forward_modelling.Forward_model._remesh_crust_lith_asth(
#         coord, np.array([105]))
#     assert np.allclose(result, np.array(
#         [0, 105, 200, 300, 400, 500, 600])) == True
#     result = forward_modelling.Forward_model._remesh_crust_lith_asth(
#         coord, np.array([5]))
#     assert np.allclose(result, np.array(
#         [0, 5, 100, 200, 300, 400, 500, 600])) == True
#     result = forward_modelling.Forward_model._remesh_crust_lith_asth(
#         coord, np.array([150]))
#     assert np.allclose(result, np.array(
#         [0, 100, 150, 200, 300, 400, 500, 600])) == True


# def test_advection():
#     f_model = instantiate_test_model()
#     a = np.arange(0, 1000, 100)
#     b = np.arange(100, 1100, 100)
#     t = np.arange(0, 100, 10)
#     r = f_model._advection(a, t, b)
#     chk = np.arange(0, 90, 10, dtype=np.float64)
#     chk = np.append(0, chk)
#     assert np.allclose(r, chk) == True
#     r = f_model._advection(b, t, a)
#     assert np.allclose(r, np.arange(10, 110, 10, dtype=np.float64) ) == True


# def test_initial_conditions():
#     f_model = instantiate_test_model()
#     f_model._setup_initial_conditions()
#     # test self._heat_production
#     with open('./tests/data/initial_HP.pickle', 'rb') as file:
#         check = pickle.load(file)
#     assert np.allclose(f_model.current_node.initial_crustal_HP, check) == True
#     #assert f_model.current_node.initial_crustal_HP == 1
#     # test self._update_initial_kLith()
#     assert round(f_model.current_node.kLith, 1) == 3.2
#     # test self._initial_temperature()
#     with open('./tests/data/t_init.pickle', 'rb') as file:
#         check = pickle.load(file)
#     assert np.allclose(f_model.current_node.Tinit, check) == True
#     # test self._initial_height_of_sealevel()
#     assert round(f_model.current_node.H0, 1) == 12768.9


# def test_update_lithosphere_depth():
#     f_model = instantiate_test_model()
#     T_arr = np.array([5, 10, 500, 1000, 1500, 2000])
#     coord = np.array([0, 100, 200, 500, 600, 1000])
#     depth_LAB = f_model._update_lithosphere_depth(T_arr, coord)
#     assert depth_LAB == 563.34
#     T_arr = np.array([5, 10, 500, 1000, 1316.7, 2000])
#     coord = np.array([0, 100, 200, 500, 600, 1000])
#     depth_LAB = f_model._update_lithosphere_depth(T_arr, coord)
#     assert depth_LAB == 600


# def test_coord_rift_scaler():
#     pass


# def test_distribute_beta_factor():
#     pass


# def test_simulate_continental():
#     pass


# def test_sediments_mean_porosity():
#     pass


# def test_sedimentation():
#     pass


# def test_compaction():
#     pass


# def test_decompaction():
#     pass


# def test_add_sediments():
#     pass


# def test_get_all_sediemnts():
#     pass


# def test_recompact_old_sediments():
#     pass
