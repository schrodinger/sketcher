import { expect, test } from '@playwright/test';
import { Sketcher } from './sketcher/wrappers/sketcher.js';
import { widgetRect } from './sketcher/wrappers/sketcher_wasm.js';

test.describe('tst_help_menu', () => {
  test('main', async ({ page }) => {
    const sk = new Sketcher(page);
    await sk.open();

    // Source: sk.help_menu('getting_started'). The legacy test relies on
    // Squish exception handling to dismiss this dialog; click its visible OK
    // button here so the next user-visible Help action can run.
    await sk.help_menu('getting_started');
    expect((await widgetRect(page, 'ok_btn')).width).toBeGreaterThan(0);
    await sk.click_button('ok_btn');

    // Source: sk.help_menu('about_sketcher'); sk.click_button('about_sketcher_ok')
    await sk.help_menu('about_sketcher');
    await sk.click_button('about_sketcher_ok');
  });
});
