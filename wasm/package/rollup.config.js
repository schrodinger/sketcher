import commonjs from '@rollup/plugin-commonjs';
import md5file from 'md5-file';
import peerDepsExternal from 'rollup-plugin-peer-deps-external';
import replace from '@rollup/plugin-replace';
import resolve from '@rollup/plugin-node-resolve';
import typescript from '@rollup/plugin-typescript';

const packageJson = require('./package.json');

export default {
  input: 'src/index.ts',
  // react/jsx-runtime is pulled in by the automatic JSX transform; keep it (and
  // react) external alongside peerDepsExternal so it is not bundled.
  external: [/^react($|\/)/],
  output: [
    {
      file: packageJson.main,
      format: 'cjs',
      sourcemap: true,
    },
    {
      file: packageJson.module,
      format: 'esm',
      sourcemap: true,
    },
  ],
  plugins: [
    // The iframe loads the sketcher shell as `static/wasm_shell.html?cache_bust=<hash>`.
    // The hash busts the shell itself; it is computed from the already-assembled
    // WASM asset that the publish flow moves into dist/ before this build runs.
    replace({
      preventAssignment: true,
      __WASM_HASH__: `${md5file.sync('dist/wasm_shell.html')}`,
    }),
    peerDepsExternal(),
    resolve({ extensions: ['.mjs', '.js', '.json', '.node', '.ts', '.tsx'] }),
    commonjs(),
    typescript({
      tsconfig: './tsconfig.json',
      declaration: true,
      declarationDir: 'dist/react',
      outDir: 'dist/react',
    }),
  ],
};
