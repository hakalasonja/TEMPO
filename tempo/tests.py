"""
Created on Wed Nov 20

This file includes all test cases for the TEMPO library to be run with pytest.

@author:
"""
import numpy as np
import pytest
from qutip import Qobj, sigmax, sigmay, sigmaz, basis, Options, expect
import math

from tempo.hamiltonian import Hamiltonian
from tempo.pulse_recipe import Pulse_recipe
from tempo.pulse import Pulse
from tempo.pulse_sequence import Pulse_sequence
from tempo.evolver import Evolver
from tempo.exceptions import *



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

class TestFullSystem:

    def test_JJ(self):
        
        # 2-level system initial state
        state_init = basis(2,1) 

        # Applied magnetic field
        Bx = 0; By = 0; Bz = 1000 # B field components (G)
        B0 = np.array([Bx, By, Bz]) # B field vector (G)

        # helper function for Hamiltonian terms
        def dotproduct(vecV, vecU):
            #
            # Dot product between vector V = (Vx,Vy,Vz) and vector U = (Ux,Uy,Uz)
            #
            return sum([Vcomp*Ucomp for Vcomp, Ucomp in zip(vecV, vecU)])


        # define function for static Hamiltonian
        def func_Zeeman(Hmats, Hpars):
            return Hpars['coeff']*Hpars['gamma']*dotproduct(Hpars['Bfield'], Hmats['S1'])

        # parameters for static Hamiltonian 
        pars_Zeeman = {'coeff': -2*np.pi, 'Bfield': B0, 'gamma': -2.8025}
        mats_Zeeman = {'S1': (sigmax(), sigmay(), sigmaz())}

        # create static Hamiltonian object
        H_Zeeman = Hamiltonian(mats_Zeeman, pars_Zeeman, func_Zeeman)

        # get the static Hamiltonian operator
        Hstat = H_Zeeman.H

        # calculate eigenenergies
        eig_energies = Hstat.eigenenergies()/2/np.pi

        # resonance frequency of transition
        frq_trans =  max(eig_energies) - min(eig_energies)
                
        # define drive amplitude 
        gamma = -2.8025 # electron gyromagnetic ratio
        Omega = 30   # on-resonant Rabi frequency
        B_amp = Omega / gamma # amplitude of magnetic field 

        # define the time-dependent pulse type of fucnctional form:
        # A*cos(omega*(t+offset)+phase)
        # function must have inputs t, args
        def func_X(t, args):
            return args['B_amp']*np.cos(2*np.pi*args['freq']*t)

        # operator to be multiplied by function output
        H_X = 2 * np.pi * gamma * sigmax()
        
        # parameter names for pulse recipe
        keys_X = ['B_amp', 'freq']

        # create the pulsetype object with the operator, parameter names, and function
        # we are not inputting numerical parameters yet; this pulsetype is a blueprint for any ACD pulse
        # values will come in later
        recipe_X = Pulse_recipe(Hamiltonian(H_X), keys_X, func_X)

        # create a pulse sequence
        pseq = Pulse_sequence(Hstat = Hstat)

        # simulation time
        tpulse = .1 

        # create array of evaluation times for output
        npts_eval = 51
        times_eval = np.linspace(0, tpulse, npts_eval)

        t0 = 0 # pulse start time

        pulseX = Pulse(recipe_X, start_time = t0, duration = tpulse, coeff_params = {'B_amp': B_amp, 'freq': frq_trans})
        pseq.add_pulse(pulseX)

        # define solver options
        opts = Options(rhs_reuse = False, nsteps = 1000000, atol=1e-9, rtol=1e-9, tidy=False)

        # initialize evolver objects (no difference at this point)
        ev = Evolver(state_init, times_eval, pseq, opts = opts)

        # Run evaluation
        result = ev.evolve(method = 'regular')
        
        expectations_Z = expect(sigmaz(), result.states)
        
        expected_result_Z =  [-1.0,                   -0.929289167227064,     -0.7305134890711701,    -0.4233591334991757,    
                              -0.06520482126415034,   0.3105095378983825,     0.6371636330303323,     0.8758275976559109,
                              0.992369285264469,      0.9692342286621923,     0.8075202090657745,     0.5373768494538111,     
                              0.18672953649056434,    -0.1880238037326829,    -0.5342753330335834,    -0.8105018963653297,    
                              -0.9679197364087414,    -0.9918544221210958,    -0.8767714784069418,    -0.6376880059957346, 
                              -0.30752019892974564,   0.0603667302309191,     0.4281886123510632,     0.7274148932986982,     
                              0.9302446845226064,     0.9999964185183823,     0.9302509040466652,     0.7274265344102241,
                              0.4282039417033618,     0.060383661411712086,   -0.3075040547088368,    -0.6376749414273759,
                              -0.8767633237162767,    -0.9918522637037576,    -0.9679239961923038,    -0.8105118338968311,    
                              -0.5342896802122861,    -0.1880404790132824,    0.18671286302881834,    0.5373625317541852, 
                              0.8075101983878971,     0.9692300538399772,     0.9923713756209621,     0.8758357842500599,     
                              0.6371767176957583,     0.3105256705768797,     -0.06518788167009171,   -0.4233437552879566,    
                              -0.7305019029266913,    -0.9292829035348242,    -0.9999999998496966 
                              ]
        
        atol = 1e-7
        rtol = 1e-7
        for i, x in enumerate(expected_result_Z):
            assert(math.isclose(x, expectations_Z[i], abs_tol=atol, rel_tol=rtol))
