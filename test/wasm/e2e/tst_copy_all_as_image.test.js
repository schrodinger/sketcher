import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';

// SKETCH-2851 certification: the top-bar More Actions path must update the
// clipboard with a structure image. This intentionally does not use the
// context menu, whose image action was already wired before the fix.
const ATOMISTIC_STRUCTURE = 'c1ccccc1O';

test.describe('SKETCH-2851 Copy All As Image', () => {
  test('top-bar menu copies an atomistic structure as a PNG', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    await test.step('Create an atomistic structure through the visible Import menu', async () => {
      await sk.import_menu('paste_in_text', ATOMISTIC_STRUCTURE);
    });

    await test.step('More Actions -> Copy All As -> Image writes a pasteable PNG', async () => {
      const image = await sk.copy_all_as_image();
      await expect(image).toMatchObject({ type: 'image/png' });
      await expect(image.size).toBeGreaterThan(100);
      await expect(image.width).toBeGreaterThan(20);
      await expect(image.height).toBeGreaterThan(20);
    });
  });
});
