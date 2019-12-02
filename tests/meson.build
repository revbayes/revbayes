run_integration_tests_sh = find_program('./run_integration_tests.sh')

test('run_integration_tests',
     run_integration_tests_sh,
     timeout: 3000,
     workdir: join_paths(meson.source_root(), 'tests'),
     args: [rb.full_path()])

