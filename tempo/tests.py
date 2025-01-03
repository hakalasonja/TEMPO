"""
Created on Wed Nov 20

This file includes all test cases for the TEMPO library to be run with pytest.

@author:
"""
import pytest
import numpy as np
import types
from collections.abc import Iterable

from qutip import sigmax, sigmay, sigmaz, basis, Options

from tempo.exceptions import TEMPO_ImproperInputException, TEMPO_NullValueException
from tempo.hamiltonian import Hamiltonian
from tempo.pulse_recipe import Pulse_recipe
from tempo.pulse import Pulse
from tempo.pulse_sequence import Pulse_sequence
from tempo.evolver import Evolver


'''
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######

                            hamiltonian.py

##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
'''

class TestHamiltonian:
    """
    Tests for the Hamiltonian class, focusing on edge cases with improper inputs.
    """

    def test_initialization_ops_none(self):
        with pytest.raises(TEMPO_ImproperInputException, match="Matrices must be in a dictionary"):
            Hamiltonian(ops=None)

    def test_initialization_ops_invalid_type(self):
        with pytest.raises(TEMPO_ImproperInputException, match="Matrices must be in a dictionary"):
            Hamiltonian(ops="invalid")

    def test_initialization_ops_qobj_without_func(self):
        # Valid case where ops is a Qobj and func is not provided
        ham = Hamiltonian(ops=sigmax())
        assert isinstance(ham, Hamiltonian)

    def test_initialization_ops_dict_without_func(self):
        ops = {'H1': sigmax()}
        with pytest.raises(TEMPO_ImproperInputException, match="Must be a function"):
            Hamiltonian(ops=ops)

    def test_initialization_func_invalid_type(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0}
        with pytest.raises(TEMPO_ImproperInputException):
            Hamiltonian(ops=ops, op_params=op_params, func="not a function")

    def test_initialization_op_params_invalid_type(self):
        ops = {'H1': sigmax()}
        func = lambda ops, params: params['a'] * ops['H1']
        with pytest.raises(TEMPO_ImproperInputException):
            Hamiltonian(ops=ops, op_params="not a dict", func=func)

    def test_initialization_op_params_contains_non_qobj(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0, 'b': 'not a float'}
        func = lambda ops, params: params['a'] * ops['H1']
        Hamiltonian(ops=ops, op_params=op_params, func=func)

    def test_initialization_correct(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0}
        func = lambda ops, params: params['a'] * ops['H1']
        ham = Hamiltonian(ops=ops, op_params=op_params, func=func)
        assert isinstance(ham, Hamiltonian)

    def test_set_ops_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException, match="Matrices must be in a dictionary"):
            ham.ops = "invalid"

    def test_set_op_params_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException):
            ham.op_params = "invalid"

    def test_set_func_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException, match="Must be a function"):
            ham.func = "invalid"

    def test_H_no_func_defined(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0}
        with pytest.raises(TEMPO_ImproperInputException):
            ham = Hamiltonian(ops=ops, op_params=op_params)
            _ = ham.H

    def test_H_func_returns_non_qobj(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0}
        func = lambda ops, params: "not a Qobj"
        ham = Hamiltonian(ops=ops, op_params=op_params, func=func)
        _ = ham.H

    def test_set_H_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException, match="The full Hamiltonian must be a QuTiP Qobj object."):
            ham.H = "invalid"

    def test_delete_ops(self):
        ham = Hamiltonian(ops=sigmax())
        del ham.ops
        with pytest.raises(AttributeError):
            _ = ham.ops

    def test_delete_op_params(self):
        ham = Hamiltonian(ops=sigmax())
        del ham.op_params
        with pytest.raises(AttributeError):
            _ = ham.op_params

    def test_delete_func(self):
        ham = Hamiltonian(ops=sigmax())
        del ham.func
        with pytest.raises(AttributeError):
            _ = ham.func

    def test_delete_H(self):
        ops = {'H1': sigmax()}
        op_params = {'a': 1.0}
        func = lambda ops, params: params['a'] * ops['H1']
        ham = Hamiltonian(ops=ops, op_params=op_params, func=func)
        ham.H  # Generate H to ensure it's set
        _ = ham.H

'''
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######

                            pulse_recipe.py

##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
'''

class TestPulseRecipe:
    """
    Tests for the Pulse_recipe class, focusing on edge cases with improper inputs.
    """

    def test_initialization_ham_none(self):
        Pulse_recipe(ham=None, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)

    def test_initialization_ham_invalid_type(self):
        with pytest.raises(TEMPO_ImproperInputException, match="Hamiltonian is not a hamiltonian object"):
            Pulse_recipe(ham="invalid", param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)

    def test_initialization_param_keys_none(self):
        ham = Hamiltonian(ops=sigmax())
        Pulse_recipe(ham=ham, param_keys=None, coeff_func=lambda t, args: args['amp'] * t)

    def test_initialization_param_keys_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        # (note this should pass because technically a string is an iterable, and we COULD encode 
        # keys one by one in a string)
        Pulse_recipe(ham=ham, param_keys="invalid", coeff_func=lambda t, args: args['amp'] * t)
        with pytest.raises(TEMPO_ImproperInputException, match="Param keys are not an iterable"):
            Pulse_recipe(ham=ham, param_keys=15, coeff_func=lambda t, args: args['amp'] * t)

    def test_initialization_coeff_func_none(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException):
            Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=None)

    def test_initialization_coeff_func_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        with pytest.raises(TEMPO_ImproperInputException, match="Coefficient function is not callable"):
            Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func="invalid")

    def test_initialization_correct(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: args['amp'] * t
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=coeff_func)
        assert isinstance(pulse_recipe, Pulse_recipe)

    def test_eval_coeff_invalid_t_type(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: args['amp'] * t
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=coeff_func)
        with pytest.raises(TypeError):
            pulse_recipe.eval_coeff(t="invalid", params={'amp': 1.0})

    def test_eval_coeff_params_not_dict(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: args['amp'] * t
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=coeff_func)
        with pytest.raises(TypeError):
            pulse_recipe.eval_coeff(t=0.0, params="invalid")

    def test_eval_coeff_params_missing_keys(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: args['amp'] * t
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp', 'freq'], coeff_func=coeff_func)
        pulse_recipe.eval_coeff(t=0.0, params={'amp': 1.0})

    def test_eval_coeff_return_type_mismatch(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: "invalid return type"
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=coeff_func)
        with pytest.raises(TEMPO_ImproperInputException, match="The type of t must match the return type of the provided coefficient function"):
            pulse_recipe.eval_coeff(t=0.0, params={'amp': 1.0})

    def test_set_ham_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        with pytest.raises(TEMPO_ImproperInputException, match="ham must be a Hamiltonian instance"):
            pulse_recipe.ham = "invalid"

    def test_set_param_keys_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        with pytest.raises(TEMPO_ImproperInputException):
            pulse_recipe.param_keys = 15

    def test_set_coeff_func_invalid_type(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        with pytest.raises(TEMPO_ImproperInputException, match="Must be a function"):
            pulse_recipe.coeff_func = "invalid"

    def test_delete_ham(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        del pulse_recipe.ham
        with pytest.raises(AttributeError):
            _ = pulse_recipe.ham

    def test_delete_param_keys(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        del pulse_recipe.param_keys
        with pytest.raises(AttributeError):
            _ = pulse_recipe.param_keys

    def test_delete_coeff_func(self):
        ham = Hamiltonian(ops=sigmax())
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp'], coeff_func=lambda t, args: args['amp'] * t)
        del pulse_recipe.coeff_func
        with pytest.raises(AttributeError):
            _ = pulse_recipe.coeff_func


'''
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######

                            pulse.py

##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
'''

class TestPulse:
    """
    Tests for the Pulse class, focusing on edge cases with improper inputs.
    """

    def test_initialization_recipe_none(self):
        with pytest.raises(TEMPO_NullValueException, match="Recipe for pulse is null"):
            Pulse(None)

    def test_initialization_recipe_invalid_type(self):
        with pytest.raises(TEMPO_ImproperInputException, match="Recipe is not a pulse_recipe.Pulse_recipe"):
            Pulse("invalid")

    def test_initialization_start_time_invalid(self):
        recipe = self._get_valid_recipe()
        with pytest.raises(TEMPO_ImproperInputException, match="Start time is not a float"):
            Pulse(recipe, start_time="invalid")

    def test_initialization_duration_invalid(self):
        recipe = self._get_valid_recipe()
        with pytest.raises(TEMPO_ImproperInputException):
            Pulse(recipe, duration="invalid")

    def test_initialization_coeff_params_invalid(self):
        recipe = self._get_valid_recipe()
        with pytest.raises(TEMPO_ImproperInputException, match="Coefficient parameters are not a dictionary"):
            Pulse(recipe, coeff_params="invalid")

    def test_initialization_correct(self):
        recipe = self._get_valid_recipe()
        coeff_params = {'param1': 1.0, 'param2': 2.0}
        pulse = Pulse(recipe, start_time=0.0, duration=1.0, coeff_params=coeff_params)
        assert isinstance(pulse, Pulse)

    def test_update_params_keys_not_iterable(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Keys are not an iterable"):
            pulse.update_params(15, [1.0])

    def test_update_params_values_not_iterable(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Values are not an iterable"):
            pulse.update_params(['key1'], 15)

    def test_update_params_keys_include_non_string(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Keys include a non string"):
            pulse.update_params(['key1', 2], [1.0, 2.0])

    def test_update_params_values_include_non_float(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Values include a non float"):
            pulse.update_params(['key1', 'key2'], [1.0, 'not a float'])

    def test_eval_coeff_t_invalid(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Time cannot be cast to float"):
            pulse.eval_coeff('invalid')

    def test_eval_coeff_args_not_dict(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Args were not provided in dictionary form"):
            pulse.eval_coeff(0.0, args='not a dict')

    def test_set_recipe_invalid_type(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="pulserecipe must be a member of the Pulserecipe class"):
            pulse.recipe = "invalid"

    def test_set_start_time_invalid(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Time could not be cast to float"):
            pulse.start_time = "invalid"

    def test_set_duration_invalid(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException):
            pulse.duration = "invalid"

    def test_set_end_time_before_start_time(self):
        pulse = self._get_valid_pulse()
        pulse.start_time = 1.0
        with pytest.raises(TEMPO_ImproperInputException, match="Start time must be before end time"):
            pulse.end_time = 0.5

    def test_set_coeff_params_invalid_type(self):
        pulse = self._get_valid_pulse()
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse coefficient parameters must be in a dictionary"):
            pulse.coeff_params = "invalid"

    def test_set_coeff_params_missing_keys(self):
        pulse = self._get_valid_pulse()
        coeff_params = {'missing_param': 1.0}
        with pytest.raises(KeyError):
            pulse.coeff_params = coeff_params

    def test_delete_recipe(self):
        pulse = self._get_valid_pulse()
        del pulse.recipe
        with pytest.raises(AttributeError):
            _ = pulse.recipe

    def test_delete_start_time(self):
        pulse = self._get_valid_pulse()
        del pulse.start_time
        with pytest.raises(AttributeError):
            _ = pulse.start_time

    def test_delete_duration(self):
        pulse = self._get_valid_pulse()
        del pulse.duration
        with pytest.raises(AttributeError):
            _ = pulse.duration

    def test_delete_end_time(self):
        pulse = self._get_valid_pulse()
        del pulse.end_time
        with pytest.raises(AttributeError):
            _ = pulse.end_time

    def test_delete_coeff_params(self):
        pulse = self._get_valid_pulse()
        del pulse.coeff_params
        with pytest.raises(AttributeError):
            _ = pulse.coeff_params

    def test_delete_ham(self):
        pulse = self._get_valid_pulse()
        del pulse.ham
        with pytest.raises(AttributeError):
            _ = pulse.ham

    # Helper methods to create valid instances
    def _get_valid_recipe(self):
        ham = Hamiltonian(ops=sigmax())
        param_keys = ['param1', 'param2']
        coeff_func = lambda t, args: args['param1'] * np.sin(args['param2'] * t)
        recipe = Pulse_recipe(ham=ham, param_keys=param_keys, coeff_func=coeff_func)
        return recipe

    def _get_valid_pulse(self):
        recipe = self._get_valid_recipe()
        coeff_params = {'param1': 1.0, 'param2': 2.0}
        pulse = Pulse(recipe=recipe, start_time=0.0, duration=1.0, coeff_params=coeff_params)
        return pulse


'''
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######

                            evolver.py

##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
'''

class TestEvolver:
    """
    Tests for the Evolver class, focusing on edge cases with improper inputs.
    """

    def test_initialization_state_init_none(self):
        Evolver(state_init=None)

    def test_initialization_state_init_invalid_type(self):
        with pytest.raises(TEMPO_ImproperInputException, match="Initial state is not a Qobj"):
            Evolver(state_init="invalid")

    def test_initialization_tlist_none(self):
        state_init = basis(2, 0)
        Evolver(state_init=state_init, tlist=None)

    def test_initialization_tlist_invalid_type(self):
        state_init = basis(2, 0)
        with pytest.raises(TEMPO_ImproperInputException):
            Evolver(state_init=state_init, tlist="invalid")

    def test_initialization_tlist_contains_non_float(self):
        state_init = basis(2, 0)
        tlist = [0.0, 'invalid', 1.0]
        with pytest.raises(TEMPO_ImproperInputException):
            Evolver(state_init=state_init, tlist=tlist)

    def test_initialization_pulse_seq_none(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse sequence is not a Pulse Sequence object"):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq="BETA")

    def test_initialization_pulse_seq_invalid_type(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse sequence is not a Pulse Sequence object"):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq="invalid")

    def test_initialization_c_ops_invalid_type(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        pulse_seq = self._get_valid_pulse_sequence()
        with pytest.raises(TEMPO_ImproperInputException):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq=pulse_seq, c_ops="invalid")

    def test_initialization_c_ops_contains_non_qobj(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        pulse_seq = self._get_valid_pulse_sequence()
        c_ops = [sigmax(), "invalid"]
        with pytest.raises(TEMPO_ImproperInputException, match="c_ops includes a non Qobj"):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq=pulse_seq, c_ops=c_ops)

    def test_initialization_e_ops_invalid_type(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        pulse_seq = self._get_valid_pulse_sequence()
        with pytest.raises(TEMPO_ImproperInputException):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq=pulse_seq, e_ops="invalid")

    def test_initialization_e_ops_contains_invalid_elements(self):
        state_init = basis(2, 0)
        tlist = [0.0, 1.0]
        pulse_seq = self._get_valid_pulse_sequence()
        e_ops = [sigmax(), 15]
        with pytest.raises(TEMPO_ImproperInputException):
            Evolver(state_init=state_init, tlist=tlist, pulse_seq=pulse_seq, e_ops=e_ops)

    def test_generate_H_pulse_seq_invalid_type(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException):
            evolver.generate_H(pulse_seq="invalid")

    def test_serial_generate_H_pulselist_invalid_type(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException):
            evolver.serial_generate_H(pulselist="invalid")

    def test_serial_generate_H_pulselist_contains_non_pulse(self):
        evolver = self._get_valid_evolver()
        pulselist = [self._get_valid_pulse(), "invalid"]
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse list includes a non Pulse object"):
            evolver.serial_generate_H(pulselist=pulselist)

    def test_serial_generate_H_safe_invalid_type(self):
        evolver = self._get_valid_evolver()
        pulselist = [self._get_valid_pulse()]
        with pytest.raises(TEMPO_ImproperInputException, match="safe is not a bool"):
            evolver.serial_generate_H(pulselist=pulselist, safe="invalid")

    def test_evolve_method_invalid(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException):
            evolver.evolve(method="invalid")

    def test_evolve_tlist_contains_non_float(self):
        evolver = self._get_valid_evolver()
        tlist = [0.0, 'invalid', 1.0]
        with pytest.raises(TEMPO_ImproperInputException, match="time list includes a non-float value"):
            evolver.evolve(tlist=tlist)

    def test_evolve_tlist_values_too_close(self):
        evolver = self._get_valid_evolver()
        tlist = [0.0, 1-(1e-9), 1.0]
        with pytest.raises(TEMPO_ImproperInputException, match="tlist values too close together; try decreasing t_rtol or changing tlist"):
            evolver.evolve(tlist=tlist, method='serial')

    def test_evolve_pulse_seq_invalid_type(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse sequence is not a Pulse Sequence object"):
            evolver.evolve(pulse_seq="invalid")

    def test_evolve_c_ops_contains_non_qobj(self):
        evolver = self._get_valid_evolver()
        c_ops = [sigmax(), "invalid"]
        with pytest.raises(TEMPO_ImproperInputException, match="c_ops includes a non Qobj"):
            evolver.evolve(c_ops=c_ops)

    def test_evolve_e_ops_contains_invalid_elements(self):
        evolver = self._get_valid_evolver()
        e_ops = [sigmax(), "invalid"]
        with pytest.raises(TEMPO_ImproperInputException):
            evolver.evolve(e_ops=e_ops)

    def test_set_pulse_seq_invalid_type(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException, match="Pulse sequence must be a Pulse_sequence instance"):
            evolver.pulse_seq = "invalid"

    def test_set_Hstat_invalid_type(self):
        evolver = self._get_valid_evolver()
        with pytest.raises(TEMPO_ImproperInputException, match="Static Hamiltonian operator must be a Hamiltonian instance or a QuTiP Qobj instance"):
            evolver.Hstat = "invalid"

    def test_delete_state_init(self):
        evolver = self._get_valid_evolver()
        del evolver.state_init
        with pytest.raises(AttributeError):
            _ = evolver.state_init

    def test_delete_tlist(self):
        evolver = self._get_valid_evolver()
        del evolver.tlist
        with pytest.raises(AttributeError):
            _ = evolver.tlist

    def test_delete_pulse_seq(self):
        evolver = self._get_valid_evolver()
        del evolver.pulse_seq
        with pytest.raises(AttributeError):
            _ = evolver.pulse_seq

    def test_delete_Hstat(self):
        evolver = self._get_valid_evolver()
        del evolver.Hstat
        with pytest.raises(AttributeError):
            _ = evolver.Hstat

    def test_delete_c_ops(self):
        evolver = self._get_valid_evolver()
        del evolver.c_ops
        with pytest.raises(AttributeError):
            _ = evolver.c_ops

    def test_delete_e_ops(self):
        evolver = self._get_valid_evolver()
        del evolver.e_ops
        with pytest.raises(AttributeError):
            _ = evolver.e_ops

    # Helper methods to create valid instances
    def _get_valid_evolver(self):
        state_init = basis(2, 0)
        tlist = np.linspace(0, 1, 10)
        pulse_seq = self._get_valid_pulse_sequence()
        evolver = Evolver(state_init=state_init, tlist=tlist, pulse_seq=pulse_seq)
        return evolver

    def _get_valid_pulse_sequence(self):
        pulse = self._get_valid_pulse()
        pulse_seq = Pulse_sequence(pulses=[pulse])
        return pulse_seq

    def _get_valid_pulse(self):
        ham = Hamiltonian(ops=sigmax())
        coeff_func = lambda t, args: args['amp'] * np.sin(args['freq'] * t)
        pulse_recipe = Pulse_recipe(ham=ham, param_keys=['amp', 'freq'], coeff_func=coeff_func)
        coeff_params = {'amp': 1.0, 'freq': 2.0}
        pulse = Pulse(recipe=pulse_recipe, start_time=0.0, duration=1.0, coeff_params=coeff_params)
        return pulse


'''
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######

                            full-use tests.

##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ###### ##### ######
'''
from qutip import *
import numpy as np
import scipy as scipy
import scipy.signal as scs
from multiprocess import Pool
import time
from scipy.optimize import curve_fit

import sys
import os

from tempo.hamiltonian import Hamiltonian
from tempo.pulse_recipe import Pulse_recipe
from tempo.evolver import Evolver
from tempo.pulse_sequence import Pulse_sequence
from tempo.pulse import Pulse
from tempo.hamfuncs import *

import math
from qutip import identity, jmat, tensor, basis, qobj
import os
import numpy as np

from tempo.exceptions import *
from collections.abc import Iterable

class TestFullSystem:
    class Qsys:
        """
        A class for representing a system of one or more spins and their spin operators.
        
        The object is constructed by passing a tuple of the dimensions of each spin's Hilbert space. For example, for a coupled system of a spin-1 and spin-1/2, use ``dimensions = (3,2)``.
        
        Also includes methods for obtaining the spin operators Sx, Sy, and Sz, as well as a list of the basis vectors for each spin.
        
        Parameters
        ----------
        dimensions : tuple
            A tuple of the dimensions of each spin's Hilbert space.
            The length of the tuple is the number of spins in the system.
            
        Attributes
        ----------
        dimensions: tuple of int
            A tuple of the dimensions of each spin's Hilbert space.
            The length of the tuple is the number of spins in the system.
        numparticles : int
            Number of spins in the system.
        stot : list of float
            A list of the spin quantum numbers of each spin. 
        Sx : list of `qutip.Qobj`
            A list of the Sx spin operators of each spin.
        Sy : list of `qutip.Qobj`
            A list of the Sy spin operators of each spin.
        Sz : list of `qutip.Qobj`
            A list of the Sz spin operators of each spin.
        """
        
        def __init__(self, dimensions: tuple):
            """
            Qsys constructor.

            Parameters
            ----------
            dimensions: tuple of int
                A tuple of the dimensions of each spin's Hilbert space.
                The length of the tuple is the number of spins in the system.
            """
            self._dimensions = dimensions

            # Input check
            if not isinstance(self.dimensions, tuple) or len(dimensions) == 0:
                raise TEMPO_ImproperInputException("Dimensions is not a tuple, or is empty")

            for c_idx, obj in enumerate(self.dimensions):
                if not isinstance(obj, int):
                    raise TEMPO_ImproperInputException(f"Dimensions includes a non-int object at index {c_idx}")
                if obj <= 0:
                    raise TEMPO_ImproperInputException("Invalid Dimension")
                

            self._numparticles = len(self._dimensions)
            
            self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
            
            self._Sx = self.oprs('x')
            self._Sy = self.oprs('y')
            self._Sz = self.oprs('z')
            
        def oprs(self, axis: str) -> Iterable[qobj.Qobj]:
            """
            Spin operators of each spin along `axis`.
            
            Parameters
            ----------
            axis : str {'x', 'y', 'z'}
                Axis along which to create spin operators; `'x'`, `'y'`, or `'z'`.
                
            Returns
            ----------
            oprs : list of `qutip.Qobj`
                List of (`qutip.Qobj`) spin operators along `axis` in the same order as spins in the `dimensions` tuple. 
            """
            if axis != "x" and axis != "y" and axis != "z":
                raise TEMPO_ImproperInputException("Axis for operator creation is not x,y,z")

            oprs = []
            idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
            
            for i in np.arange(self._numparticles):
                idlist[i] = jmat(self._stot[i], axis)
                oprs.append(tensor(idlist))
                idlist[i] = identity(self._dimensions[i])
            
            #oprs.append(oprs[0]+oprs[1])
            
            return oprs
        
        def basisstates(self, particlenum: int) -> Iterable[qobj.Qobj]:
            """
            Basis state vectors of spin at index `particlenum` in the `dimensions` tuple. 
            
            Parameters
            ----------
            particlenum : int
                Index of spin for which to generate basis vectors.
            
            Returns
            ----------
            basisstates : list of `qutip.Qobj`
                List of (`qutip.Qobj`) basis state vectors. 
                Length of list is the dimension of the spin's Hilbert space.
            """
            if not castable(particlenum, int):
                raise TEMPO_ImproperInputException("particlenum is not an integer")
            if particlenum >= len(self._dimensions) or particlenum < -1 * len(self._dimensions):
                raise TEMPO_ImproperInputException(f"Spin index {particlenum} was given a value outside the range -1 * len(dimensions) --> len(dimensions) - 1.")

            particlenum = int(particlenum)

            basisstates = []
            
            idlist = [identity(self._dimensions[i]) for i in np.arange(self._numparticles)]
            
            for i in np.arange(self._dimensions[particlenum]):
                idlist[particlenum] = basis(self._dimensions[particlenum], i)*basis(self._dimensions[particlenum], i).dag()
                basisstates.append(tensor(idlist))
                idlist[particlenum] = identity(self._dimensions[particlenum])
            
            return basisstates
        
        @property
        def dimensions(self):
            return self._dimensions
        
        @dimensions.setter
        def dimensions(self, dimensions):
            self._dimensions = dimensions
            self._numparticles = len(self._dimensions)
            self._stot = [(self._dimensions[i]-1)/2 for i in np.arange(self._numparticles)]
            self._Sx = self.oprs('x')
            self._Sy = self.oprs('y')
            self._Sz = self.oprs('z')

        
        @dimensions.deleter 
        def dimensions(self):
            del self._dimensions
        
        @property
        def numparticles(self):
            return self._numparticles
        
        @numparticles.setter
        def numparticles(self, num):
            self._numparticles = num
        
        @numparticles.deleter 
        def numparticles(self):
            del self._numparticles
        
        @property
        def stot(self):
            return self._stot
        
        @stot.setter
        def stot(self, ls):
            self._stot = ls
        
        @stot.deleter 
        def stot(self):
            del self._stot
        
        @property
        def Sx(self):
            return self._Sx
        
        @Sx.setter
        def Sx(self, Sx):   
            if type(Sx) == qobj.Qobj:
                self._Sx = Sx
        
        @Sx.deleter 
        def Sx(self):
            del self._Sx
        
        @property
        def Sy(self):
            return self._Sy
        
        @Sy.setter
        def Sy(self, Sy): 
            if type(Sy) == qobj.Qobj:
                self._Sy = Sy
        
        @Sy.deleter 
        def Sy(self):
            del self._Sy
        
        @property
        def Sz(self):
            return self._Sz
        
        @Sz.setter       
        def Sz(self, Sz): 
            if type(Sz) == qobj.Qobj:
                self._Sz = Sz

        @Sz.deleter   
        def Sz(self):
            del self._Sz

    def test_ramsey_1(self):
        def ret_transition_frq(Ham, idx_mi0, idx_mif, idx_ms0=1, idx_msf=0):
            #
            # Return frequency difference (MHz) between two eigenstates 
            #

            eigenvals, eigenstates = Ham.eigenstates()
            state0 = tensor(basis(qs.dimensions[0], idx_ms0), basis(qs.dimensions[1], idx_mi0))
            statef = tensor(basis(qs.dimensions[0], idx_msf), basis(qs.dimensions[1], idx_mif))

            idx_e0 = np.argmax([np.sum(np.abs(eigenstates[i].full())*state0.full()) for i in np.arange(len(eigenvals))])
            idx_ef = np.argmax([np.sum(np.abs(eigenstates[i].full())*statef.full()) for i in np.arange(len(eigenvals))])

            frq_trans = (eigenvals[idx_ef] - eigenvals[idx_e0])/2/np.pi 
            return frq_trans #MHz
        
        # initialize the NV system: (3, 2) for the dimensions of the electron's and nucleus' Hilbert spaces
        # creates a coupled system of the two spins
        qs = TestFullSystem.Qsys((3,2))
        state_init = tensor(basis(3,1), basis(2, 0)) # ms = 0 (electron), mi = +1/2 state (nucleus)

        #
        # Applied magnetic field
        #
        Bx = 0; By = 0; Bz = 300 # B field components (G)
        B0_init = np.array([Bx, By, Bz]) # B field vector (G)

        # If desired, rotate magnetic field around y-axis, by angle theta
        theta = 0 #degrees
        theta *= np.pi/180
        Ry = np.array([[np.cos(theta), 0, np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]]) # rotation opr
        B0 = np.dot(Ry, B0_init) 
        # create Hamiltonian objects to store these terms
        HZFS = ZFS(qs)
        HZeeNV = Zeeman(qs, B0)
        HZeeNuc = Zeeman(qs, B0, nuc = True)
        HHF = Hyperfine(qs)

        # get the total static Hamiltonian operator
        Hstat = HZFS.H + HZeeNV.H + HZeeNuc.H + HHF.H
        
        frq_trans = ret_transition_frq(Hstat, 0, 0, idx_ms0=1, idx_msf=0) # ms=0 <=> ms=+1 transition
        assert(frq_trans == 3710.991165214304)

        # Number corresponds approximately to rabi freq in MHz (exact if no detuning)

        gammaNV = -2.8025

        BacMHz = 30   # AC B-field frequency, ms=+1 transition, units MHz

        B_amp = BacMHz / gammaNV # amplitude of B-field 

        detunings = [0.5, 1, 2] # ms = +1 detuning (MHz)


        # From rabi frequency & detuning, predict generalized Rabi frequency & determine length of pi pulse    
        frqs_rabi = [((BacMHz)**2 + abs(det)**2)**0.5 for det in detunings]
            
        # Define duration of pi/2 pulse from predicted Rabi frequency    
        t_pipulse = [1/frq/2 for frq in frqs_rabi]

        #
        # Array of times for free evolution
        # Keep only values greater than duration of pi/2 pulse to avoid overlapping pulses
        #
        tfp_max = 2 # us
        arr_timesfp = np.linspace(0, tfp_max, 41)
        arr_timesfp = arr_timesfp[arr_timesfp > t_pipulse[2]/2]


        # define the time-dependent pulse type
        # in this case AC drive, A*cos(omega*(t+offset)+phase)
        # function must have inputs t, args
        def ACD(t, args):
            return args['B_amp']*np.cos(2*np.pi*args['freq']*(t + args['offset']) + args['phase'])

        H_Bac = 2 * np.pi * gammaNV * qs.Sx[0] * np.sqrt(2)
        ACDkeys = ['B_amp', 'freq', 'offset', 'phase']

        # create the pulsetype object with the operator, parameter names, and function
        # we are not inputting numerical parameters yet; this pulsetype is a blueprint for any ACD pulse
        # values will come in later
        ACDpulsetype = Pulse_recipe(Hamiltonian(H_Bac), ACDkeys, ACD)

        ls_pulsehandlers_Ramsey = []
        tlists = []

        for i in np.arange(len(detunings)):

            # for each time duration of free evolution we generated, run the sequence
            for t_fp in arr_timesfp:
            
                # pulse sequence handler object
                duration_1 = t_pipulse[i]/2
                pulse_handler_Ramsey = Pulse_sequence(Hstat = Hstat)
                tlists.append(np.linspace(-duration_1/2, t_fp+duration_1, 100))

                
                # Loop through two Ramsey pulses
                for tally_p in np.arange(2):
                
                    # define pulse timings
                    starttime_1 = tally_p*t_fp - duration_1/2
                
                    Bamp = B_amp
                    frq = frq_trans+detunings[i]
                    offset = 0
                    phase = 0
                    
                    pulse_ramsey_90 = Pulse(ACDpulsetype, start_time = starttime_1, duration = duration_1, coeff_params = {'B_amp': Bamp, 'freq': frq, 'offset': offset, 'phase': phase})
                    pulse_handler_Ramsey.add_pulse(pulse_ramsey_90)
                        
                # add pulse sequence handlers to list    
                ls_pulsehandlers_Ramsey.append(pulse_handler_Ramsey)

        # define how to execute pulse
        # this creates an Evolver object and then calls its evolve() function, which returns the state of the system at the timestamps specified in tlist
        def exec_pulse_Ramsey(state_init, tlist, pulsehandler):
            opts = Options(rhs_reuse = False, nsteps = 1000000, atol=1e-9, rtol=1e-9, tidy=False)
            return Evolver(state_init, tlist, pulsehandler, opts = opts).evolve(method = 'serial')

        # execute pulse sequences in parallel
        sim_starttime = time.time()

        n = len(ls_pulsehandlers_Ramsey)

        # tlist has to be sufficiently large for qutip's mesolve to be able to use it
        # inputs are initial state, tlist, and pulse sequence, as outlined in the args for exec_pulse_Ramsey
        inputs = zip([state_init]*n, tlists, ls_pulsehandlers_Ramsey)

        with Pool(processes = n) as pool:
            output_map = pool.starmap(exec_pulse_Ramsey, inputs)
            pool.close()
            pool.join()

        # output_map = list(output_map)

        # if you are on Mac, and want to see a progress bar, you can uncomment below
        # output_map = parallel_map(exec_pulse_Ramsey, np.arange(len(ls_pulsehandlers_Ramsey)), 
        #                             task_args=(state_init, tlists, ls_pulsehandlers_Ramsey), progress_bar=True)

        sim_endtime = time.time()

        # the only purpose of the code in this cell is to separate the output of the parallel map into the 
        # three detunings. The first 40 entries in output_parallelmap are for detuning = 0.5 MHz, etc. 

        n_detunings = len(detunings) # number of separate time sweeps, one for each detuning
        len_timesweep = int(len(output_map)/n_detunings) # length of one time-sweep array's results

        ls_states_timesweep = [[]]*n_detunings # create a 2D array to store each of the sweeps' results

        for i in np.arange(n_detunings):
            ls_states_timesweep[i] = [elem.states[-1] for elem in output_map[i*len_timesweep:(i+1)*len_timesweep]] 
            # By taking the very last element of each pulse, we get the state of the system at the end of 
            # each pulse

        popdata = [[]]*n_detunings
        exp_op = qs.basisstates(0)[1]

        expected_res = [np.array([0.00753587, 0.05681842, 0.15033208, 0.28058748, 0.42797377,
            0.58270007, 0.73386795, 0.85885509, 0.94749185, 0.99479686,
            0.99245684, 0.94308892, 0.84981552, 0.71951621, 0.57180232,
            0.41734384, 0.26631612, 0.1410456 , 0.0524556 , 0.00522694]),
            np.array([0.02817681, 0.21253892, 0.51165658, 0.80131285, 0.97830604,
            0.97181222, 0.78743444, 0.48845009, 0.19858431, 0.02171102,
            0.02821734, 0.2125158 , 0.51149282, 0.80146421, 0.97823406,
            0.97178005, 0.78744481, 0.48862725, 0.19843558, 0.0217811 ]),
            np.array([0.1061628 , 0.67031253, 0.99972584, 0.63679219, 0.0870632 ,
            0.10621763, 0.67019323, 0.99972153, 0.63683917, 0.0870906 ,
            0.10627034, 0.67007257, 0.9997172 , 0.6368941 , 0.08711588,
            0.10632   , 0.66995135, 0.99971277, 0.63695461, 0.08713781])]

        for i in np.arange(n_detunings):
            popdata[i] = expect(exp_op, ls_states_timesweep[i])
            c = 0
            for j in range(len(popdata[i]), 2):
                assert(math.isclose(popdata[i][j], expected_res[i][c], abs_tol=1e-7, rel_tol=1e-7))
                c += 1
    
    def test_ramsey_2(self):
        def ret_transition_frq(Ham, idx_mi0, idx_mif, idx_ms0=1, idx_msf=0):
            #
            # Return frequency difference (MHz) between two eigenstates 
            #

            eigenvals, eigenstates = Ham.eigenstates()
            state0 = tensor(basis(qs.dimensions[0], idx_ms0), basis(qs.dimensions[1], idx_mi0))
            statef = tensor(basis(qs.dimensions[0], idx_msf), basis(qs.dimensions[1], idx_mif))

            idx_e0 = np.argmax([np.sum(np.abs(eigenstates[i].full())*state0.full()) for i in np.arange(len(eigenvals))])
            idx_ef = np.argmax([np.sum(np.abs(eigenstates[i].full())*statef.full()) for i in np.arange(len(eigenvals))])

            frq_trans = (eigenvals[idx_ef] - eigenvals[idx_e0])/2/np.pi 
            return frq_trans #MHz
        
        # initialize the NV system: (3, 2) for the dimensions of the electron's and nucleus' Hilbert spaces
        # creates a coupled system of the two spins
        qs = TestFullSystem.Qsys((3,2))
        state_init = tensor(basis(3,1), basis(2, 0)) # ms = 0 (electron), mi = +1/2 state (nucleus)

        #
        # Applied magnetic field
        #
        Bx = 0; By = 0; Bz = 300 # B field components (G)
        B0_init = np.array([Bx, By, Bz]) # B field vector (G)

        # If desired, rotate magnetic field around y-axis, by angle theta
        theta = 0 #degrees
        theta *= np.pi/180
        Ry = np.array([[np.cos(theta), 0, np.sin(theta)],[0,1,0],[-np.sin(theta), 0, np.cos(theta)]]) # rotation opr
        B0 = np.dot(Ry, B0_init) 
        # create Hamiltonian objects to store these terms
        HZFS = ZFS(qs)
        HZeeNV = Zeeman(qs, B0)
        HZeeNuc = Zeeman(qs, B0, nuc = True)
        HHF = Hyperfine(qs)

        # get the total static Hamiltonian operator
        Hstat = HZFS.H + HZeeNV.H + HZeeNuc.H + HHF.H
        
        frq_trans = ret_transition_frq(Hstat, 0, 0, idx_ms0=1, idx_msf=0) # ms=0 <=> ms=+1 transition
        exp_op = qs.basisstates(0)[1]

        
        # Number corresponds approximately to rabi freq in MHz (exact if no detuning)

        gammaNV = -2.8025

        BacMHz = 30   # AC B-field frequency, ms=+1 transition, units MHz

        B_amp = BacMHz / gammaNV # amplitude of B-field 

        detunings = [0.5, 1, 2] # ms = +1 detuning (MHz)


        # define the time-dependent pulse type
        # in this case AC drive, A*cos(omega*(t+offset)+phase)
        # function must have inputs t, args
        def ACD(t, args):
            return args['B_amp']*np.cos(2*np.pi*args['freq']*(t + args['offset']) + args['phase'])

        # define how to execute pulse
        # this creates an Evolver object and then calls its evolve() function, which returns the state of the system at the timestamps specified in tlist
        def exec_pulse_Ramsey(state_init, tlist, pulsehandler):
            opts = Options(rhs_reuse = False, nsteps = 1000000, atol=1e-9, rtol=1e-9, tidy=False)
            return Evolver(state_init, tlist, pulsehandler, opts = opts).evolve(method = 'serial')

        H_Bac = 2 * np.pi * gammaNV * qs.Sx[0] * np.sqrt(2)
        ACDkeys = ['B_amp', 'freq', 'offset', 'phase']

        # create the pulsetype object with the operator, parameter names, and function
        # we are not inputting numerical parameters yet; this pulsetype is a blueprint for any ACD pulse
        # values will come in later
        ACDpulsetype = Pulse_recipe(Hamiltonian(H_Bac), ACDkeys, ACD)

        det = 1
        frq_rabi = ((BacMHz)**2 + abs(det)**2)**0.5
        t_pipulse = 1/frq_rabi/2
        #print(t_pipulse) # units us

        # try different pulse durations: pi/4, pi/3, pi/2, 3pi/2
        t_pulses = [t_pipulse/4, t_pipulse/3, t_pipulse/2, 3*t_pipulse/4]

        #
        # Array of times for free evolution
        # Keep only values greater than duration of pi/4 pulse to avoid overlapping pulses
        #
        tfp_max = 2 # us
        arr_timesfp = np.linspace(0, tfp_max, 41)
        arr_timesfp = arr_timesfp[arr_timesfp > t_pulses[3]]

        ls_pulsehandlers_Ramsey = []
        tlists = []

        for i in np.arange(len(t_pulses)):

            # for each time duration of free evolution we generated, run the sequence
            for t_fp in arr_timesfp:
            
                # pulse sequence handler object
                duration_1 = t_pulses[i]
                pulse_handler_Ramsey = Pulse_sequence(Hstat = Hstat)
                tlists.append(np.linspace(-duration_1/2, t_fp+duration_1, 100))

                
                # Loop through two Ramsey pulses
                for tally_p in np.arange(2):
                
                    # define pulse timings
                    starttime_1 = tally_p*t_fp - duration_1/2
                
                    Bamp = B_amp
                    frq = frq_trans+det
                    phase = 0
                    offset = 0

                    # create pulse object and add to pulse sequence handler
                    pulse_ramsey_90 = Pulse(ACDpulsetype, start_time = starttime_1, duration = duration_1, coeff_params = {'B_amp': Bamp, 'freq': frq, 'offset': offset, 'phase': phase})
                    pulse_handler_Ramsey.add_pulse(pulse_ramsey_90)
                        
                # add pulse sequence handlers to list    
                ls_pulsehandlers_Ramsey.append(pulse_handler_Ramsey)

        sim_starttime = time.time()


        n = len(ls_pulsehandlers_Ramsey)
        inputs = zip([state_init]*n, tlists, ls_pulsehandlers_Ramsey)

        with Pool(processes = n) as pool:
            output_map = pool.starmap(exec_pulse_Ramsey, inputs)
            pool.close()
            pool.join()

        output_map = list(output_map)

        sim_endtime = time.time()
        print('Time taken (s)', round(sim_endtime - sim_starttime,3))

        n_pulsedurations = len(t_pulses) # number of separate time sweeps, one for each pulse duration
        len_timesweep = int(len(output_map)/n_pulsedurations) # length of one time-sweep array's results

        ls_states_timesweep = [[]]*n_pulsedurations # create a 2D array to store each of the sweeps' results

        for i in np.arange(n_pulsedurations):
            ls_states_timesweep[i] = [elem.states[-1] for elem in output_map[i*len_timesweep:(i+1)*len_timesweep]]

        popdata = [[]]*n_pulsedurations
        colors = ['firebrick', 'gold', 'forestgreen', 'blue']

        expected_res = [
            np.array([0.51116362, 0.54810724, 0.60407261, 0.67286966, 0.74899093,
       0.8263921 , 0.89708559, 0.9528112 , 0.98796277, 0.99999674,
       0.98792814, 0.95226783, 0.89635978, 0.82635388, 0.74946904,
       0.67283902, 0.60346033, 0.54786135, 0.51153689, 0.49866914,
       0.51137728, 0.54828321, 0.60434595, 0.67320852, 0.74922286,
       0.82643966, 0.89705308, 0.95281246, 0.98798532, 0.99999639,
       0.98794374, 0.95235736, 0.8964987 , 0.82648208, 0.74959272,
       0.67302124, 0.60375195, 0.548264  , 0.51198013, 0.49903325]),
            np.array([0.27212193, 0.32399821, 0.4060142 , 0.51061952, 0.62738235,
       0.74360858, 0.84698853, 0.92839329, 0.98132388, 0.99999386,
       0.9809245 , 0.92653098, 0.84426302, 0.74193479, 0.62680273,
       0.50930086, 0.40357756, 0.32213621, 0.27170368, 0.25479499,
       0.2725701 , 0.32442503, 0.40629058, 0.51070463, 0.62736857,
       0.74365074, 0.84710759, 0.92845848, 0.9812848 , 0.99998289,
       0.98101991, 0.92659465, 0.84420041, 0.74194021, 0.62710308,
       0.50977937, 0.40391333, 0.32224134, 0.27179567, 0.25507423]),
            np.array([2.81768086e-02, 1.00969217e-01, 2.12538921e-01, 3.54469069e-01,
       5.11656576e-01, 6.65465259e-01, 8.01312849e-01, 9.09071427e-01,
       9.78306044e-01, 9.99849303e-01, 9.71812218e-01, 8.98946594e-01,
       7.87434440e-01, 6.45639740e-01, 4.88450092e-01, 3.34500287e-01,
       1.98584306e-01, 9.08893791e-02, 2.17110212e-02, 1.67439584e-04,
       2.82173360e-02, 1.01084100e-01, 2.12515805e-01, 3.54236673e-01,
       5.11492821e-01, 6.65562186e-01, 8.01464206e-01, 9.09056907e-01,
       9.78234063e-01, 9.99852244e-01, 9.71780054e-01, 8.98831260e-01,
       7.87444813e-01, 6.45866074e-01, 4.88627247e-01, 3.34420943e-01,
       1.98435581e-01, 9.08971432e-02, 2.17810981e-02, 1.66991685e-04]),
            np.array([0.51988683, 0.55882503, 0.6192722 , 0.69266564, 0.77005096,
       0.84571749, 0.91308609, 0.96394199, 0.99308258, 0.99913758,
       0.98096026, 0.93865659, 0.87814315, 0.80809936, 0.73294319,
       0.65537159, 0.5844311 , 0.53315331, 0.50726436, 0.50356107,
       0.51992064, 0.55872569, 0.61901263, 0.69249781, 0.7699872 ,
       0.84561187, 0.91298242, 0.96394295, 0.99311652, 0.99913115,
       0.98101505, 0.93875349, 0.87801088, 0.80772727, 0.73276425,
       0.65558629, 0.58458074, 0.53279418, 0.50668375, 0.50329453])
        ]

        for i in np.arange(n_pulsedurations):
            popdata[i] = expect(exp_op, ls_states_timesweep[i])
            for j in range(len(popdata[i])):
                assert(math.isclose(popdata[i][j], expected_res[i][j], abs_tol=1e-7, rel_tol=1e-7))
